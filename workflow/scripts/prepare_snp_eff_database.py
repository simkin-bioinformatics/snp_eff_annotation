import subprocess
snp_eff_folder=snakemake.input.snp_eff_folder
input_genome_fasta=snakemake.input.genome_fasta
input_genome_gff=snakemake.input.genome_gff
database=snakemake.params.database
description=snakemake.params.description
output_genome_fasta=snakemake.output.genome_fasta
output_genome_gff=snakemake.output.genome_gff
temp_config=snakemake.output.temp_config
input_config=snakemake.input.snp_eff_folder+'/snpEff.config'

subprocess.call(f'mkdir {snp_eff_folder}/data', shell=True)
subprocess.call(f'mkdir {snp_eff_folder}/data/{database}', shell=True)
subprocess.call(f'cp {input_genome_fasta} {output_genome_fasta}', shell=True)
subprocess.call(f'cp {input_genome_gff} {output_genome_gff}', shell=True)

output_file=open(temp_config, 'w')
for line in open(input_config):
	output_file.write(line)
	if 'data.dir = ./data/' in line:
		output_file.write(f'\n{database}.genome : {description}\n\n')
output_file.close()

subprocess.call(f'cp {temp_config} {input_config}', shell=True)