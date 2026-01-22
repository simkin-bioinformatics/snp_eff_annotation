'''
parses the annotations of snpEff to output an AA table file for downstream use.
The current version does not handle targeted vs. nontargeted mutations. Basic
algorithm is as follows:
1. reads through a vcf file and stores sample columns
2. upon reaching a snpEff line, iterates first through all samples, then all 15
fields of snpEff. Because a single snpEff VCF line can contain multiple
mutations, any field that has remainder 0 in modulo 15 is treated as the
beginning of a new mutation. All fields are then stored in a nested dictionary,
keyed first by sample and then by mutation. format and clair3 output are also
captured
3. Search the gff file to look up gene names associated with gene IDs and use
snpEff annotations to fill in the 'header' portion of each mutation. Parse out
the clair3 'AD' field values into cov, ref, and alt, where ref is the first 'AD'
value, alt is the second, and coverage is the ref+alt (I don't use DP for
coverage because DP is always slightly bigger than summed AD values - probably
due to reads that clair3 considers 'sequencing error').
4. convert snpEff header fields and counts to AA tables.

The 6 AA table fields are:
Gene ID
Gene
Mutation Name
ExonicFunc
AA Change
Targeted

By comparison, the 15 fields snpEff outputs are:
'Allele 
 Annotation 
 Annotation_Impact 
 Gene_Name 
 Gene_ID 
 Feature_Type 
 Feature_ID 
 Transcript_BioType 
 Rank 
 HGVS.c 
 HGVS.p 
 cDNA.pos / cDNA.length 
 CDS.pos / CDS.length 
 AA.pos / AA.length 
 Distance 

The mapping of snpEff onto AA table is:
Gene_ID->Gene ID
Annotation->ExonicFunc,
gff lookup of Name associated with Gene_ID->Gene
HGVS.p->AA Change
Gene+AA Change->Mutation Name
Targeted - currently hardcoded as 'unspecified'
'''

#input_vcf='annotated_variants.vcf'
input_vcf=snakemake.input.snp_eff_vcf
input_gff=snakemake.input.genome_gff
targets_tsv=snakemake.input.targets_tsv
coverage_AA=snakemake.output.coverage_AA
reference_AA=snakemake.output.reference_AA
alternate_AA=snakemake.output.alternate_AA
output_paths=[coverage_AA, reference_AA, alternate_AA]

def make_targets_dict(targets_tsv, gene_mappings):
	'''
	converts a targets.tsv file into a dictionary, such that the chromosome,
	reference allele, and alternate allele can be easily compared against
	entries from a vcf file. The bottom layer is a header that can be added to
	an output AA table. A later function will examine the VCF file lines and add
	sequencing depths.
	'''
	#[f'{ref+alt},{ref},{alt}', Gene_ID, gene_name, Mutation_Name, ExonicFunc, AA_change, targeted]
	targets_dict, targeted_mutations={}, []
	for line_number, line in enumerate(open(targets_tsv)):
		if line_number>0:
			line=line.strip().split('\t')
			chrom, pos, ID, ref, alt=line[0:5]
			mut=line[7]
			gene_ID=line[11]
			gene_name=gene_mappings[gene_ID]
			mutation_name=gene_name+'-'+mut
			targets_dict.setdefault(chrom, {})
			targets_dict[chrom].setdefault(pos, {})
			targets_dict[chrom][pos].setdefault(ref, {})
			targets_dict[chrom][pos][ref][alt]=[gene_ID, gene_name, mutation_name, 'missense_mutation', mut, 'Yes']
			targeted_mutations.append(mutation_name)
	return targets_dict, targeted_mutations

def extract_counts(labels, values):
	if len(labels) != len(values):
		values = "./."+":."*(len(labels) - 1)
	labels=labels.split(':')
	values=values.split(':')
	for label_number, label in enumerate(labels):
		if label=='AD':
			depths=values[label_number].split(',')
			break
	if len(depths)==2:
		ref, alt=depths
		ref=convert_count(ref)
		alt=convert_count(alt)
	elif len(depths)==1:
		ref=depths[0]
		ref=convert_count(ref)
		alt=0
	elif len(depths)>2:
		ref, alt=depths[0], depths[1]
		ref=convert_count(ref)
		alt=convert_count(alt)
	else:
		print('weird depths', labels, values)
	return ref, alt

def check_unannotated_targeted(vcf_line, targets_dict, depth_dict, samples, mutation_number):
	'''
	This function only applies to vcf entries that have no snpEff annotation.
	Checks a vcf_line to see if it's "targeted" using targets_dict and extracts
	reference and alternate counts for each sample if it is targeted. Outputs
	the resulting depths to depth_dict
	'''
	chrom, pos, ID, ref_allele, alt_allele=vcf_line[:5]
	if chrom in targets_dict and pos in targets_dict[chrom] and ref_allele in targets_dict[chrom][pos]:
		for alt_allele in targets_dict[chrom][pos][ref_allele]:
			header_list=targets_dict[chrom][pos][ref_allele][alt_allele]
			for sample_number, sample in enumerate(samples):
				labels, sample_counts=vcf_line[8], vcf_line[9+sample_number]
				ref, alt=extract_counts(labels, sample_counts)
				depth_dict.setdefault(sample, {})
				depth_dict[sample][mutation_number]=[f'{ref+alt},{ref},{alt}']+header_list
			mutation_number+=1
	return depth_dict, mutation_number

def grab_mutation_names(sample_dict):
	'''
	grabs all the mutation names from the first sample of a sample dictionary
	'''
	first_sample=list(sample_dict.keys())[0]
	return [sample_dict[first_sample][mutation][3] for mutation in sample_dict[first_sample]]

def make_zeroes_dict(samples, missing_muts):
	'''
	creates a dictionary assigning 0 coverage, ref, and alt counts to any
	mutations that were in the targeted mutations but missing from both the
	snpEff dictionary and the dictionary of targeted mutations that had vcf file
	coverage.
	'''
	zeroes_dict, mut={},1
	for chrom in targets_dict:
		for pos in targets_dict[chrom]:
			for ref in targets_dict[chrom][pos]:
				for alt in targets_dict[chrom][pos][ref]:
					if targets_dict[chrom][pos][ref][alt][2] in missing_muts:
						header=targets_dict[chrom][pos][ref][alt]
						for sample in samples:
							zeroes_dict.setdefault(sample, {})
							zeroes_dict[sample][mut]=['0,0,0']+header
						mut+=1
	return zeroes_dict

def remove_redundant(redundant_muts, eff_dict, covered_dict):
	'''
	removes any mutations in the covered_dict (vcf calls corresponding to
	targeted SNPs) that are also annotated already by snpEff - this is necessary
	because some targeted SNPs are tri-allelic, with the same position
	associated with two different mutations, only one of which might be observed
	and therefore annotated by snpEff
	'''
	new_covered={}
	first_sample=list(covered_dict.keys())[0]
	good_muts=[]
	for mut in covered_dict[first_sample]:
		if covered_dict[first_sample][mut][3] not in redundant_muts:
			for sample in covered_dict:
				new_covered.setdefault(sample, {})
				new_covered[sample][mut]=covered_dict[sample][mut]
	return new_covered

def special_sort(unsorted_list):
	'''
	sorts mutations first by gene name and then by numbered amino acid position
	'''
	sorting_list=[]
	for entry in unsorted_list:
		mutation, number=entry
		gene='-'.join(mutation.split('-')[:-1])
		#aa_pos=int(mutation.split('-')[-1][3:-3])
		numbers=mutation.split('-')[-1][3:]
		if not numbers[0].isdigit():
			print('original mutation was', mutation)
		aa_pos=''
		for char in numbers:
			if char.isdigit():
				aa_pos+=char
			else:
				break
		if aa_pos:
			aa_pos=int(aa_pos)
		sorting_list.append([gene, aa_pos, entry])
	return [item[-1] for item in sorted(sorting_list)]


def merge_dicts(eff_dict, nonzero_dict, zero_dict):
	'''
	merges together the snp_eff dictionary, the nonzero dictionary, and the zero
	dictionary. Also creates a sorted list of mutations for easy lookup in the
	final AA table
	'''
	first_sample=list(eff_dict.keys())[0]
	max_eff=max(eff_dict[first_sample].keys())
	if first_sample in nonzero_dict:
		max_nonzero=max(nonzero_dict[first_sample].keys())
	sorting_list=[]
	for mut in eff_dict[first_sample]:
		sorting_list.append((eff_dict[first_sample][mut][3], mut))
	merged_dict=eff_dict
	for sample in nonzero_dict:
		for mut1 in nonzero_dict[sample]:
			nonzero_new_mut=max_eff+mut1+1
			sorting_list.append((nonzero_dict[sample][mut1][3], nonzero_new_mut))
			merged_dict[sample][nonzero_new_mut]=nonzero_dict[sample][mut1]
		for mut2 in zero_dict[sample]:
			zero_new_mut=max_eff+max_nonzero+mut2+2
			sorting_list.append((zero_dict[sample][mut2][3], zero_new_mut))
			merged_dict[sample][zero_new_mut]=zero_dict[sample][mut2]
	sorting_list=list(set(sorting_list))
	sorted_list=special_sort(sorting_list)
	#print('sorted is', sorted_list)
	#for sample in merged_dict:
	#	print(sample, 'merged is', merged_dict[sample].keys())
	return merged_dict, sorted_list


def fill_in_missing_values(protein_dict, targets_dict, targeted_depth_dict, targeted_mutations):
	'''
	Examines all targeted_mutations from targets_dict. Any of these targeted
	mutations that don't have a snpEff annotation (are not in protein_dict) and
	that don't appear as a reference-only state in the VCF file (are not in
	targeted_depth_dict) are assumed to be of 0 coverage.
	TODO:
	Currently, anything that has a snpEff value and is targeted gets counted
	twice. Need a check to remove these 'targeted' counts and trust snpEff.
	'''
	snp_eff_muts=set(grab_mutation_names(protein_dict))
	nonzero_muts=set(grab_mutation_names(targeted_depth_dict))
	targeted_muts=set(targeted_mutations)
	missing_muts=targeted_muts-(nonzero_muts|snp_eff_muts)
	redundant_muts=nonzero_muts&snp_eff_muts
	samples=protein_dict.keys()
	zeroes_dict=make_zeroes_dict(samples, missing_muts)
	new_nonzero=remove_redundant(redundant_muts, protein_dict, targeted_depth_dict)
	merged_dict, sorting_list=merge_dicts(protein_dict, new_nonzero, zeroes_dict)
	return merged_dict, sorting_list

def parse_vcf_file(input_vcf, targets_dict):
	'''
	searches vcf file for snpEff annotations and retrieves the mutation
	information and depths associated with the snpeff annotations
	'''
	ann_dict, target_depth_dict={},{}
	mutation_number, parsed_counter=1,0
	if input_vcf.endswith('.gz'):
		import gzip
		file_handle=gzip.open(input_vcf, mode='rt')
	else:
		file_handle=open(input_vcf)
	for line_number, line in enumerate(file_handle):
		line=line.strip().split('\t')
		if line[0].startswith('#CHROM'):
			samples=line[9:]
		if not line[0].startswith('#') and len(line)>7:
			if line[7].split(';')[-1].startswith('ANN'):
				unparsed_snpeff=line[7].split('|')
				for column_number, column in enumerate(unparsed_snpeff):
					if column_number%15==0:
						header=[]
						parsed_counter+=1
					header.append(column)
					for sample_number, sample in enumerate(samples):
						ann_dict.setdefault(sample, {})
						if column_number%15==14:
							ann_dict[sample][parsed_counter]=[line[8]+';'+line[9+sample_number]]+header
			target_depth_dict, mutation_number=check_unannotated_targeted(line, targets_dict, target_depth_dict, samples, mutation_number)
	return ann_dict, target_depth_dict

def convert_count(count):
	if count!='.':
		count=int(count)
	else:
		count=0
	return count

def parse_annotations(ann_dict, gene_mappings, targeted_mutations):
	'''
	extracts only columns of interest from the annotation dictionary
	'''
	protein_dict={}
	for sample in ann_dict:
		protein_dict.setdefault(sample, {})
		for mut in ann_dict[sample]:
			columns=ann_dict[sample][mut]
			if len(columns)==16 and columns[11]: #item 11 is only populated if the mutation is protein-coding
				vcf_fields=columns[0]
				Gene_ID=columns[5]
				gene_name=gene_mappings[columns[4]]
				AA_change=columns[11][2:]
				mutation_name=gene_name+'-'+AA_change
				ExonicFunc=columns[2]
				targeted='No'
				if mutation_name in targeted_mutations:
					targeted='Yes'
				labels, values=vcf_fields.split(';')
				ref, alt=extract_counts(labels, values)
				protein_dict[sample][mut]=[f'{ref+alt},{ref},{alt}', Gene_ID, gene_name, mutation_name, ExonicFunc, AA_change, targeted]
	return protein_dict


def grab_gene_mappings(input_gff):
	'''
	looks up gene names associated with gene IDs. Here is an example line from gff file to be parsed:
	ID=PF3D7_1343700;Name=Kelch13;description=kelch protein K13;ebi_biotype=protein_coding
	'''
	gene_mappings={}
	for line in open(input_gff):
		line=line.strip().split('\t')[-1].split(';')
		if line[0].startswith('ID=') and line[1].startswith('Name='):
			ID=line[0][3:]
			name=line[1][5:]
			gene_mappings[ID]=name
		elif line[0].startswith('ID='):
			ID=line[0][3:]
			name=ID
			gene_mappings[ID]=name
	return gene_mappings

def get_targeted_status(input_tsv, snpeff_dict):
	'''
	takes a targets.tsv file as input and looks up whether each mutation in a
	snpeff_dictionary is a 'targeted' mutation of interest or not.
	'''

def format_header(output_file, merged_dict, muts):
	first_sample=list(protein_dict.keys())[0]
	header_names=['Gene ID', 'Gene', 'Mutation Name', 'ExonicFunc', 'AA Change', 'Targeted']
	for row_number, name in enumerate(header_names):
		output_line=[name]
		for mut in muts:
			output_line.append(merged_dict[first_sample][mut][row_number+1])
		output_file.write(','.join(output_line)+'\n')

def output_tables(merged_dict, sorting_list, output_paths):
	muts=[entry[1] for entry in sorting_list]
	for file_number, output_path in enumerate(output_paths):
		output_file=open(output_path, 'w')
		format_header(output_file, merged_dict, muts)
		for sample in merged_dict:
			output_line=[sample]
			for mut in muts:
				output_line.append(merged_dict[sample][mut][0].split(',')[file_number])
			output_file.write(','.join(output_line)+'\n')

gene_mappings=grab_gene_mappings(input_gff)
targets_dict, targeted_mutations=make_targets_dict(targets_tsv, gene_mappings)
ann_dict, targeted_depth_dict=parse_vcf_file(input_vcf, targets_dict)
protein_dict=parse_annotations(ann_dict, gene_mappings, targeted_mutations)
merged_dict, sorting_list=fill_in_missing_values(protein_dict, targets_dict, targeted_depth_dict, targeted_mutations)
output_tables(merged_dict, sorting_list, output_paths)