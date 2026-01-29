'''
Like the original snpEff_parser except that this is what came out when I fed the
original into gemini and asked it to refactor the code to be less convolutedand
to use 'DP' information when 'AD' is not available. It is much cleaner but it
doesn't quite parse the output AA tables correctly (missing exonicFunc and
"targeted" status rows in the header).
'''

import csv
import gzip
import sys


def get_depths(format_string, sample_string, alt_index=1):
    """
    Extracts Ref and Alt counts with fallback logic.
    Priority: 1. AD (Allele Depth) | 2. DP (Total Depth) fallback.
    """
    labels = format_string.split(":")
    values = sample_string.split(":")
    data = dict(zip(labels, values))

    ref, alt = 0, 0
    # Priority 1: AD field (specific allele depths)
    if "AD" in data and data["AD"] != ".":
        depths = data["AD"].split(",")
        try:
            ref = int(depths[0])
            # Map targeted alt to the correct VCF alt index (Edge Case 2)
            if len(depths) > alt_index and alt_index > 0:
                alt = int(depths[alt_index])
        except (ValueError, IndexError):
            pass
    # Priority 2: DP field fallback (Edge Case 3)
    elif "DP" in data and data["DP"] != ".":
        try:
            ref = int(data["DP"])
            alt = 0
        except ValueError:
            pass

    return ref, alt


def parse_ann_field(info_field):
    """
    Parses snpEff ANN field from the INFO column.
    Format: Allele|Annotation|Impact|Gene_Name|Gene_ID|...
    """
    if "ANN=" not in info_field:
        return []

    ann_part = [f for f in info_field.split(";") if f.startswith("ANN=")][0][4:]
    parsed_effects = []
    for effect in ann_part.split(","):
        fields = effect.split("|")
        if len(fields) >= 11:
            parsed_effects.append(
                {
                    "alt_allele": fields[0],
                    "gene_id": fields[4],
                    "aa_change": fields[10][2:]
                    if fields[10].startswith("p.")
                    else fields[10],
                    "exonic_func": fields[1],
                }
            )
    return parsed_effects


def load_gene_mappings(gff_path):
    """Parses GFF to map Gene IDs to common Names."""
    mapping = {}
    with open(gff_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            attrs = {
                a.split("=")[0]: a.split("=")[1]
                for a in parts[-1].split(";")
                if "=" in a
            }
            gid = attrs.get("ID")
            if gid:
                mapping[gid] = attrs.get("Name", gid)
    return mapping


def load_targets(targets_path, gene_map):
    """Loads targets_tsv into a lookup dict keyed by (chrom, pos)."""
    targets = {}
    with open(targets_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            chrom, pos = row["CHROM"], row["POS"]
            ref, alt = row["REF"], row["ALT"]
            gene_id = row["gene_id"]
            mut_aa = row["aminoacid_change"]  # e.g., N51I

            gene_name = gene_map.get(gene_id, gene_id)
            key = (chrom, pos)
            if key not in targets:
                targets[key] = []

            targets[key].append(
                {
                    "ref": ref,
                    "alt": alt,
                    "gene_id": gene_id,
                    "gene_name": gene_name,
                    "aa_change": mut_aa,
                    "mut_name": f"{gene_name}-{mut_aa}",
                }
            )
    return targets


def main():
    # File paths from Snakemake or CLI
    vcf_in = snakemake.input.snp_eff_vcf
    gff_in = snakemake.input.genome_gff
    targets_in = snakemake.input.targets_tsv

    outputs = {
        "cov": open(snakemake.output.coverage_AA, "w"),
        "ref": open(snakemake.output.reference_AA, "w"),
        "alt": open(snakemake.output.alternate_AA, "w"),
    }

    gene_map = load_gene_mappings(gff_in)
    targets = load_targets(targets_in, gene_map)
    print("targets is", targets)
    # Pre-collect all unique mutation names for the header
    all_target_muts = []
    for pos_targets in targets.values():
        for t in pos_targets:
            all_target_muts.append(t)
    # Sort logically: Gene Name then AA position
    all_target_muts.sort(
        key=lambda x: (
            x["gene_name"],
            int("".join(filter(str.isdigit, x["aa_change"])) or 0),
        )
    )

    # Write Headers
    header_rows = [
        ["Gene ID"] + [t["gene_id"] for t in all_target_muts],
        ["Gene"] + [t["gene_name"] for t in all_target_muts],
        ["Mutation Name"] + [t["mut_name"] for t in all_target_muts],
        ["AA Change"] + [t["aa_change"] for t in all_target_muts],
    ]
    for out in outputs.values():
        for row in header_rows:
            out.write(",".join(row) + "\n")

    # Data structure to hold counts per sample
    # Using a nested dict to store results as we stream the VCF
    sample_data = {}

    opener = gzip.open if vcf_in.endswith(".gz") else open
    with opener(vcf_in, "rt") as vcf:
        samples = []
        for line in vcf:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                for s in samples:
                    sample_data[s] = {}
                continue

            cols = line.strip().split("\t")
            chrom, pos, vcf_ref, vcf_alts = (
                cols[0],
                cols[1],
                cols[3],
                cols[4].split(","),
            )

            if (chrom, pos) in targets:
                ann_data = parse_ann_field(cols[7])
                format_str = cols[8]

                for s_idx, s_name in enumerate(samples):
                    sample_col = cols[9 + s_idx]

                    for target in targets[(chrom, pos)]:
                        # Determine Alt Index (Edge Case 2)
                        vcf_alt_idx = (
                            vcf_alts.index(target["alt"]) + 1
                            if target["alt"] in vcf_alts
                            else -1
                        )

                        # Match with snpEff (The Brain - Edge Case 1)
                        match = next(
                            (
                                a
                                for a in ann_data
                                if a["aa_change"] == target["aa_change"]
                            ),
                            None,
                        )

                        r_count, a_count = 0, 0
                        if match and vcf_alt_idx != -1:
                            r_count, a_count = get_depths(
                                format_str, sample_col, vcf_alt_idx
                            )
                        elif vcf_alt_idx != -1:
                            # Found at site but not annotated as this AA change (e.g., synonymous)
                            r_count, a_count = get_depths(
                                format_str, sample_col, vcf_alt_idx
                            )
                        else:
                            # Not in ALT list: get REF depth (Alt=0)
                            r_count, _ = get_depths(format_str, sample_col, 0)
                            a_count = 0

                        sample_data[s_name][target["mut_name"]] = (
                            r_count + a_count,
                            r_count,
                            a_count,
                        )

    # Final Output Writing
    for s_name in samples:
        row_cov, row_ref, row_alt = [s_name], [s_name], [s_name]
        for t in all_target_muts:
            counts = sample_data[s_name].get(
                t["mut_name"], (0, 0, 0)
            )  # Default to 0 if missing (Edge Case 1)
            row_cov.append(str(counts[0]))
            row_ref.append(str(counts[1]))
            row_alt.append(str(counts[2]))

        outputs["cov"].write(",".join(row_cov) + "\n")
        outputs["ref"].write(",".join(row_ref) + "\n")
        outputs["alt"].write(",".join(row_alt) + "\n")

    for out in outputs.values():
        out.close()


if __name__ == "__main__":
    main()
