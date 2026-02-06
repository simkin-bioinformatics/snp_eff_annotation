'''
V4: Adds calculation for 'Other' allele depths to capture multiallelic 
contributions and non-called variants at targeted sites.
'''

import csv
import gzip
import sys

def get_depths(format_string, sample_string, alt_index=1):
    """
    Extracts Ref, Alt, and 'Other' counts with dual-strategy fallback.
    Strategy 1: Sum of non-target alternates in AD field.
    Strategy 2: DP - (Ref + Alt) to capture uncalled alleles/noise.
    """
    labels = format_string.split(":")
    values = sample_string.split(":")
    data = dict(zip(labels, values))
    ref, alt, other = 0, 0, 0
    total_dp = 0
    # Process Allele Depth (AD)
    if "AD" in data and data["AD"] != ".":
        depths_str = data["AD"].split(",")
        try:
            depths = [int(d) if d != "." else 0 for d in depths_str]
            ref = depths[0]
            if len(depths) > alt_index and alt_index > 0:
                alt = depths[alt_index]
                # Other is the sum of all alternate alleles minus the one we are tracking
                other = sum(depths[1:]) - alt
            else:
                # If we are looking at Ref only or index is out of bounds for AD
                other = sum(depths[1:])
        except (ValueError, IndexError):
            pass
    elif "DP" in data and data["DP"] != ".":
            try:
                # alt and other remain 0 as per the initial assignment
                ref = int(data["DP"])
            except ValueError:
                pass
    return ref, alt, other

def parse_ann_field(info_field, gene_map):
    """
    Takes a single VCF annotation line (aka variant call) as input.
    Parses the snpEff ANN field from the INFO column. Returns a list of
    dictionaries, one dictionary for each annotation associated with the
    variant.
    """
    if "ANN=" not in info_field:
        return []
    ann_part = [f for f in info_field.split(";") if f.startswith("ANN=")][0][4:]
    parsed_effects = []
    for effect in ann_part.split(","):
        fields = effect.split("|")
        if len(fields) >= 11:
            if fields[10].startswith("p."):
                aa_change=fields[10][2:]
            else:
                aa_change=fields[10]
            parsed_effects.append(
                {"alt_allele": fields[0], "exonic_func": fields[1],
                "gene_id": fields[4], "gene_name": gene_map.get(fields[3], fields[3]),
                "aa_change": aa_change
                })
            
    #print('after parsing, SNPEff list is', parsed_effects)
    return parsed_effects

def load_gene_mappings(gff_path):
    """
    looks up gene names associated with each gene ID in the gff file. If there
    is no gene name associated with the gene ID, then the gene name is the gene
    ID
    """
    gene_mappings = {}
    for line in open(gff_path):
        line = line.strip().split("\t")[-1].split(";")
        if line[0].startswith("ID=") and line[1].startswith("Name="):
            ID = line[0][3:]
            name = line[1][5:]
            gene_mappings[ID] = name
        elif line[0].startswith("ID="):
            ID = line[0][3:]
            name = ID
            if ID not in gene_mappings:
                gene_mappings[ID] = name
    return gene_mappings

def load_targets(targets_path, gene_map):
    '''
    parses the targets.tsv file - every targeted mutation becomes a dictionary.
    A given genomic location may have multiple mutations and therefore multiple
    dictionaries associated with it.
    '''
    targets = {}
    with open(targets_path, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            chrom, pos = row["CHROM"], row["POS"]
            ref, alt = row["REF"], row["ALT"]
            gene_id = row["gene_id"]
            mut_aa = row["aminoacid_change"]

            gene_name = gene_map.get(gene_id, gene_id)
            # row_values are the extracted values from a line of targets.tsv
            row_values = {"ref": ref, "alt": alt, "gene_id": gene_id,
                "gene_name": gene_name, "aa_change": mut_aa,
                "mut_name": f"{gene_name}-{mut_aa}", "targeted": "Yes",
            }
            key = (chrom, pos)
            targets.setdefault(key, []).append(row_values)
    return targets

def reformat_targets(targets):
    '''
    reformats the targets.tsv dictionary produced from the function above so
    that instead of being keyed by genomic position, it's keyed by the mutation
    name. load_targets and reformat_targets could probably be combined into a
    single function, except that the output of load_targets is used elsewhere in
    the program.
    '''
    reformatted_targets = {}
    for pos_targets in targets.values():
        for column in pos_targets:
            reformatted_targets[column["mut_name"]] = {
                "gene_id": column["gene_id"],
                "gene_name": column["gene_name"],
                "exonic_func": "missense_variant",  # Default for targets
                "aa_change": column["aa_change"],
                "targeted": "Yes",
            }
    return reformatted_targets

def get_counts(sample_data, samples, vcf_alt_idx, cols, format_str, mut_name):
    # Extract counts for all samples
    for s_idx, s_name in enumerate(samples):
        r_count, a_count, o_count = get_depths(
            format_str, cols[9 + s_idx], vcf_alt_idx
        )
        # Store as 4-tuple: (Cov, Ref, Alt, Other)
        sample_data[s_name][mut_name] = (
            r_count + a_count + o_count, 
            r_count,                     
            a_count,                     
            o_count                      
        )
    return sample_data

def parse_annotations(ann_data, vcf_alts, observed_mutations, samples, sample_data, cols, format_str):
    """
    Iterates through snpEff annotations. If synonymous nucleotide changes result 
    in the same amino acid change, their alternate counts are summed.
    """
    for ann in ann_data:
        if not ann["aa_change"]:
            continue
        mut_name = f"{ann['gene_name']}-{ann['aa_change']}"
        try:
            vcf_alt_idx = vcf_alts.index(ann["alt_allele"]) + 1
        except ValueError:
            continue
        # 1. Update metadata tracking
        if mut_name not in observed_mutations:
            observed_mutations[mut_name] = {
                "gene_id": ann["gene_id"],
                "gene_name": ann["gene_name"],
                "exonic_func": ann["exonic_func"],
                "aa_change": ann["aa_change"],
                "targeted": "No",
            }
        else:
            observed_mutations[mut_name]["exonic_func"] = ann["exonic_func"]
        # 2. Extract counts for this specific alternate allele
        for s_idx, s_name in enumerate(samples):
            ref, alt, other = get_depths(format_str, cols[9 + s_idx], vcf_alt_idx)
            if mut_name in sample_data[s_name]:
                # SUMMING LOGIC: Add new Alt counts to existing Alt counts
                # Subtract these counts from 'Other' since they are no longer 'Other'
                curr_cov, curr_ref, curr_alt, curr_other = sample_data[s_name][mut_name]
                new_alt = curr_alt + alt
                new_other = max(0, curr_other - alt) # Ensure we don't go below zero
                sample_data[s_name][mut_name] = (curr_cov, curr_ref, new_alt, new_other)
            else:
                # FIRST ENCOUNTER: Standard 4-tuple assignment
                sample_data[s_name][mut_name] = (ref + alt + other, ref, alt, other)
    return observed_mutations, sample_data

def assign_targeted_depths(targets, targeted_key, ann_data, samples, sample_data, format_str, cols):
    """
    Handles targeted mutations not found in the VCF annotations.
    Uses get_depths to prioritize AD values.
    """
    if targeted_key in targets:
        # Create a set of mutations actually found at this site in the VCF
        annotated_muts = set([f"{a['gene_name']}-{a['aa_change']}" for a in ann_data])
        
        for target in targets[targeted_key]:
            m_name = target["mut_name"]
            # If a targeted mutation is NOT present in the VCF annotations
            if m_name not in annotated_muts:
                for s_idx, s_name in enumerate(samples):
                    if m_name not in sample_data[s_name]:
                        # Strategy: Call get_depths with an invalid alt_index (-1).
                        # This forces all alternate AD counts into the 'other' category.
                        r_count, _, o_count = get_depths(format_str, cols[9 + s_idx], alt_index=-1)
                        # Coverage is the sum of high-quality AD reads (Ref + Other)
                        # The specific target Alt is 0 because it wasn't called.
                        total_cov = r_count + o_count
                        sample_data[s_name][m_name] = (total_cov, r_count, 0, o_count)
    return sample_data

def special_sort(unsorted_list):
    import re
    sorting_list = []
    for mutation in unsorted_list:
        gene = "-".join(mutation.split("-")[:-1])
        mut_part = mutation.split("-")[-1]
        number_match = re.search(r'\d+', mut_part)
        aa_pos = int(number_match.group()) if number_match else 0
        sorting_list.append((gene, aa_pos, mutation))
    return [item[2] for item in sorted(sorting_list)]

def write_output(cov_path, ref_path, alt_path, other_path, observed_mutations, sorted_mut_names, samples, sample_data):
    output_files = {
        "cov": open(cov_path, "w"), 
        "ref": open(ref_path, "w"),
        "alt": open(alt_path, "w"),
        "other": open(other_path, "w")
    }
    
    headers = [
        ["Gene ID"] + [observed_mutations[mut]["gene_id"] for mut in sorted_mut_names],
        ["Gene"] + [observed_mutations[mut]["gene_name"] for mut in sorted_mut_names],
        ["Mutation Name"] + [mut for mut in sorted_mut_names],
        ["ExonicFunc"] + [observed_mutations[mut]["exonic_func"] for mut in sorted_mut_names],
        ["AA Change"] + [observed_mutations[mut]["aa_change"] for mut in sorted_mut_names],
        ["Targeted"] + [observed_mutations[mut]["targeted"] for mut in sorted_mut_names],
    ]

    key_map = {"cov": 0, "ref": 1, "alt": 2, "other": 3}

    for key, f in output_files.items():
        for h_row in headers:
            f.write(",".join(h_row) + "\n")
        for s_name in samples:
            row = [s_name]
            for m_name in sorted_mut_names:
                counts = sample_data[s_name].get(m_name, (0, 0, 0, 0))
                row.append(str(counts[key_map[key]]))
            f.write(",".join(row) + "\n")
        f.close()


def parse_vcf_file(vcf_in, observed_mutations, gene_map, targets):
    sample_data = {}
    opener = gzip.open if vcf_in.endswith(".gz") else open
    with opener(vcf_in, "rt") as vcf:
        samples = []
        for line in vcf:
            if line.startswith("##"): continue
            if line.startswith("#CHROM"):
                samples = line.strip().split("\t")[9:]
                for s in samples: sample_data[s] = {}
                continue
            cols = line.strip().split("\t")
            chrom, pos, vcf_alts = cols[0], cols[1], cols[4].split(",")
            ann_data = parse_ann_field(cols[7], gene_map)
            format_str = cols[8]
            observed_mutations, sample_data = parse_annotations(ann_data, vcf_alts, observed_mutations, samples, sample_data, cols, format_str)
            sample_data = assign_targeted_depths(targets, (chrom, pos), ann_data, samples, sample_data, format_str, cols)
    return observed_mutations, samples, sample_data

def main():
    # Variable names matching snpEff_parser_gemini_v3.py style
    vcf_in = snakemake.input.snp_eff_vcf
    gff_in = snakemake.input.genome_gff
    targets_in = snakemake.input.targets_tsv
    
    cov_path = snakemake.output.coverage_AA
    ref_path = snakemake.output.reference_AA
    alt_path = snakemake.output.alternate_AA
    other_path = snakemake.output.other_AA # New 4th table output

    gene_map = load_gene_mappings(gff_in)
    targets = load_targets(targets_in, gene_map)
    observed_mutations = reformat_targets(targets)

    observed_mutations, samples, sample_data = parse_vcf_file(vcf_in, observed_mutations, gene_map, targets)    
    sorted_mut_names = special_sort(list(observed_mutations.keys()))

    write_output(cov_path, ref_path, alt_path, other_path, observed_mutations, sorted_mut_names, samples, sample_data)

if __name__ == "__main__":
    main()