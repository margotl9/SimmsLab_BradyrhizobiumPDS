import pandas as pd

# Load sample information and validate
configfile: "config/config.yaml"

samples_df = pd.read_csv("config/samples.tsv", sep="\t")

# Load config paths
resources = config["resources"]
input_dir = config["input_dir"]
results = config['results']


# -----------------------------------------------------
# Phage Defense System Identification rules
# -----------------------------------------------------
# -----------------------------------------------------
# 00 Gene Annotation
# -----------------------------------------------------
rule prodigal:
	input:
		input_dir + "{location_site_plant}.fasta"
	output:
		features_gff = resources + '01_gene_protein_prediction/{location_site_plant}.gff',
		genome_faa = resources + '01_gene_protein_prediction/{location_site_plant}.faa',
	conda:
		"../envs/prodigal.yaml"
	shell:
		"""
		prodigal -i {input} -o {output.features_gff} -f gff -a {output.genome_faa}
		"""


# -----------------------------------------------------
# 01 DefenseFinder
# -----------------------------------------------------
rule defenseFinder:
	input:
		input_dir + "{location_site_plant}.fasta"
	params:
		out_dir_sp = results + '01_PHAGE_DEFENSE_SYS/{location_site_plant}/DefenseFinder/'
	output:
		results + '01_PHAGE_DEFENSE_SYS/{location_site_plant}/DefenseFinder/{location_site_plant}_defense_finder_systems.tsv'
	conda:
		"../envs/defense_finder.yaml"
	shell:
		"""
		defense-finder update
		defense-finder run {input} -o {params.out_dir_sp}
		"""

rule analyze_PDS:
	input:
		results + '01_PHAGE_DEFENSE_SYS/{location_site_plant}/DefenseFinder/{location_site_plant}_defense_finder_systems.tsv'
	output:
		results + '01_PHAGE_DEFENSE_SYS/defense_finder_analysis/{location_site_plant}_defense_finder_analysis.csv'
	params:
		df_path = resources + "df-List_Systems.csv"
	threads:
		1
	conda:
		"../envs/jupyter.yaml"
	notebook:
		"../notebooks/defense_finder_analysis.py.ipynb"


# -----------------------------------------------------
# 02 PADLOC
# -----------------------------------------------------
rule padloc_db:
	output:
		config["home_dir"] + ".snakemake/conda/e9ed3c60150a9296743ff09672bfdad2_/data/hmm/padlocdb.hmm"
	conda:
		"../envs/padloc.yaml"
	shell:
		"""
		padloc --db-update
		"""

rule padloc:
	input:
		features_gff = resources + '01_gene_protein_prediction/{location_site_plant}.gff',
		genome_faa = resources + '01_gene_protein_prediction/{location_site_plant}.faa',
		padloc_db = config["home_dir"] + ".snakemake/conda/e9ed3c60150a9296743ff09672bfdad2_/data/hmm/padlocdb.hmm"
	output:
		results + '01_PHAGE_DEFENSE_SYS/{location_site_plant}/PADLOC/{location_site_plant}_padloc.csv'
	params:
		out_dir = results + '01_PHAGE_DEFENSE_SYS/{location_site_plant}/PADLOC'
	conda:
		"../envs/padloc.yaml"
	shell:
		"""
		padloc --faa {input.genome_faa} --gff {input.features_gff} --outdir {params.out_dir} --fix-prodigal 
		"""