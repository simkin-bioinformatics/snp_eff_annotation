'''
like v2 but uses regular expressions as part of special_sort.
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
	if "AD" in data and data["AD"] != ".":
		depths = data["AD"].split(",")
		try:
			ref = int(depths[0])
			if len(depths) > alt_index and alt_index > 0:
				alt = int(depths[alt_index])
		except (ValueError, IndexError):
			pass
	elif "DP" in data and data["DP"] != ".":
		try:
			ref = int(data["DP"])
			alt = 0
		except ValueError:
			pass
	return ref, alt


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
		r_count, a_count = get_depths(
			format_str, cols[9 + s_idx], vcf_alt_idx
		)
		sample_data[s_name][mut_name] = (
			r_count + a_count,
			r_count,
			a_count,
		)
	return sample_data


def parse_annotations(ann_data, vcf_alts, observed_mutations, samples, sample_data, cols, format_str):
	'''
	iterates through all the snpEff mutations associated with the variant. If
	the mutation has an amino acid associated with it, find the position of the
	alternate allele in the AD field and get the count associated with the
	corresponding alternate allele position.
	'''
	for ann in ann_data:
		if not ann["aa_change"]:
			continue
		mut_name = f"{ann['gene_name']}-{ann['aa_change']}"
		vcf_alt_idx = vcf_alts.index(ann["alt_allele"]) + 1
		# Update metadata if not already a target (prioritize snpEff info)
		if mut_name not in observed_mutations:
			observed_mutations[mut_name] = {
				"gene_id": ann["gene_id"],
				"gene_name": ann["gene_name"],
				"exonic_func": ann["exonic_func"],
				"aa_change": ann["aa_change"],
				"targeted": "No",
			}
		else:
			# If it's a target, update the exonic function from snpEff
			observed_mutations[mut_name]["exonic_func"] = ann["exonic_func"]
		sample_data=get_counts(sample_data, samples, vcf_alt_idx, cols, format_str, mut_name)
	return observed_mutations, sample_data

def assign_targeted_depths(targets, targeted_key, ann_data, samples, sample_data, format_str, cols):
	if targeted_key in targets:
		annotated_muts=set([f"{annotation['gene_name']}-{annotation['aa_change']}" for annotation in ann_data])
		for target in targets[targeted_key]:
			if target['mut_name'] not in annotated_muts:
				for s_idx, s_name in enumerate(samples):
					# Site is target but this specific Alt is missing: use Ref depth
					if target["mut_name"] not in sample_data[s_name]:
						r_count, _ = get_depths(format_str, cols[9 + s_idx], 0)
						sample_data[s_name][target["mut_name"]] = (r_count,
							r_count, 0,)
	return sample_data


def special_sort(unsorted_list):
	"""
	Sorts mutations first by gene name and then by the first numerical 
	amino acid position found in the mutation name.
	"""
	import re
	sorting_list = []
	for mutation in unsorted_list:
		# Split by '-' and join all but the last part to get the gene name
		# Handles gene names that might contain hyphens
		gene = "-".join(mutation.split("-")[:-1])
		
		# Get the mutation part (e.g., 'N51I')
		mut_part = mutation.split("-")[-1]
		
		# Use regex to find the first sequence of digits in the mutation part
		# This is more robust than hardcoded slicing
		number_match = re.search(r'\d+', mut_part)
		
		if number_match:
			aa_pos = int(number_match.group())
		else:
			# Fallback for mutations with no numbers (e.g., 'Gene-Unknown')
			aa_pos = 0
			
		sorting_list.append((gene, aa_pos, mutation))
	
	# Sort by Gene (index 0) then Position (index 1)
	# Returns only the original mutation string (index 2)
	return [item[2] for item in sorted(sorting_list)]

def write_output(cov_path, ref_path, alt_path, observed_mutations, sorted_mut_names, samples, sample_data):
	output_files = {"cov": open(cov_path, "w"), "ref": open(ref_path, "w"),
	"alt": open(alt_path, "w")}
	# 6-row Header Construction
	headers = [
		["Gene ID"] + [observed_mutations[mut]["gene_id"] for mut in sorted_mut_names],
		["Gene"] + [observed_mutations[mut]["gene_name"] for mut in sorted_mut_names],
		["Mutation Name"] + [mut for mut in sorted_mut_names],
		["ExonicFunc"]
		+ [observed_mutations[mut]["exonic_func"] for mut in sorted_mut_names],
		["AA Change"] + [observed_mutations[mut]["aa_change"] for mut in sorted_mut_names],
		["Targeted"] + [observed_mutations[mut]["targeted"] for mut in sorted_mut_names],
	]
	for key, f in output_files.items():
		for h_row in headers:
			f.write(",".join(h_row) + "\n")
		for s_name in samples:
			row = [s_name]
			for m_name in sorted_mut_names:
				counts = sample_data[s_name].get(m_name, (0, 0, 0))
				idx = 0 if key == "cov" else (1 if key == "ref" else 2)
				row.append(str(counts[idx]))
			f.write(",".join(row) + "\n")
		f.close()

def parse_vcf_file(vcf_in, observed_mutations, gene_map, targets):
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
				cols[0], cols[1], cols[3], cols[4].split(","),
				)
			ann_data = parse_ann_field(cols[7], gene_map)
			format_str = cols[8]

			# Step 1: Process every annotation in the VCF (found mutations)
			observed_mutations, sample_data=parse_annotations(ann_data, vcf_alts, observed_mutations, samples, sample_data, cols, format_str)

			# Step 2: Handle targets at this position not found as alternates in this VCF line
			targeted_key=(chrom, pos)
			sample_data=assign_targeted_depths(targets, targeted_key, ann_data, samples, sample_data, format_str, cols)
	return observed_mutations, samples, sample_data


def main():
	vcf_in = snakemake.input.snp_eff_vcf
	gff_in = snakemake.input.genome_gff
	targets_in = snakemake.input.targets_tsv
	cov_path=snakemake.output.coverage_AA
	ref_path=snakemake.output.reference_AA
	alt_path=snakemake.output.alternate_AA

	gene_map = load_gene_mappings(gff_in)
	targets = load_targets(targets_in, gene_map)

	# Track all observed mutations: {mut_name: {metadata}}
	observed_mutations = reformat_targets(targets)

	observed_mutations, samples, sample_data=parse_vcf_file(vcf_in, observed_mutations, gene_map, targets)	

	# Sorting logic for columns
	sorted_mut_names = special_sort(list(observed_mutations.keys()))

	# Output writing
	write_output(cov_path, ref_path, alt_path, observed_mutations, sorted_mut_names, samples, sample_data)

if __name__ == "__main__":
	main()
