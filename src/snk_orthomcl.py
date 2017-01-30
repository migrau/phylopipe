#########################################################################################
## OrthoMCL Snakemake script                                                            #
# It runs all the steps of orthoMCL automatically, parallelizing blast.                 #
# Steps are based on http://orthomcl.org/common/downloads/software/v2.0/UserGuide.txt   #
# Input files: CDS proteins (obtanied from transdecoder for example).                   #
# Output file: groups.txt, containing all the ortholog clusters.                        #
#                                                                                       #
## Considerations                                                                       #
# - An empty mysql database is required.                                                #
# - FASTA files should be renamed using only four letters: XXXX.fasta                   #
# - orthomcl.config file is also required (attached).                                   #
# - Blast step is parallelized, by default 50 threads.                                  #
#                                                                                       #
## Possible error:                                                                      #
# 1. For rule adjustFASTA, the user has to indicate the position of the unique ID.      #
# By default, it is 2, change if necessary.                                             #
#                                                                                       #
## Usage                                                                                #
# $ snakemake --snakefile snakeorthomcl.py --config fastadir=../data/ -j 24 (-np)       #
#                                                                                       #
#########################################################################################

import subprocess,sys
from os.path import join

# Globals ---------------------------------------------------------------------

# Full path to a folder that holds all of your FASTA files.
FASTA_DIR = config["fastadir"]

# A Snakemake regular expression matching the FASTA files.
SAMPLES, = glob_wildcards(join(FASTA_DIR, '{sample,[^/]+}.fasta'))

# Patterns using the 'sample' wildcard.
PATTERN = '{sample}.fasta'

blastJobs=50

SAMPLESX=[]
for i in range(0,blastJobs):
	SAMPLESX.append("goodProteins."+"%02d"%i)

#SAMPLESX = ["file.0","file.1"]
# Functions -------------------------------------------------------------------

# Rules -----------------------------------------------------------------------

rule all:
	input:
		"groups.txt"

rule installSchema:
	input:
		'orthomcl.config'
	output:
		'install_schema.log'
	shell:"""
		orthomclInstallSchema {input} {output}
	"""

rule adjustFASTA:
	input:
		files=join(FASTA_DIR, PATTERN),
		db='install_schema.log'
	output:
		compliantFiles='compliantFasta/{sample}.fasta'
	shell:"""
		name=`echo -n {input.files} | tail -c10 | head -c4`;
		orthomclAdjustFasta $name {input.files} 2;
		mv $name.fasta compliantFasta/;
	"""

rule filterFASTA:
	input: 
		compliantFiles=expand('compliantFasta/{sample}.fasta', sample=SAMPLES)
	output:
		goodPro="goodProteins.fasta",
		poorPro="poorProteins.fasta"
	shell:"""
		orthomclFilterFasta compliantFasta/ 10 20
	"""

rule splitFASTA:
	input:
		"goodProteins.fasta"
	output:
		trimmed=expand('blast/{sample}.fasta', sample=SAMPLESX),
		datab='blast/prot_db.psq'
	shell:"""
		pyfasta split -n {blastJobs} {input};
		makeblastdb -in {input} -dbtype prot -out blast/prot_db;
		mv goodProteins.*.fasta blast/;
	"""
#makeblastdb -in {input} -dbtype prot -out test_data/prot_db

FASTA_DIR2 = 'blast/'
SAMPLES2 = glob_wildcards(join(FASTA_DIR2, '{sample2,[^/]+}.fasta'))
PATTERN2 = '{sample2}.fasta'

rule blast:
 	input:
 		fastas=join(FASTA_DIR2, PATTERN2),
 		datab='blast/prot_db.psq'
 	output:
 		'blastout/{sample2}.tsv'
 	shell:"""
 		blastp -db blast/prot_db -query {input.fastas} -outfmt 6 -out {output} -num_threads 8
 	"""

rule catblast:
 	input:
 		expand("blastout/{sample}.tsv", sample=SAMPLESX)
 	output:
 		'blastout/merged.tsv'
 	shell:"""
 		cat {input} > {output}
 	"""

rule BlastParser:
	input:
		blastfile='blastout/merged.tsv',
		compliantFiles=expand('compliantFasta/{sample}.fasta', sample=SAMPLES)
	output:
		"similarSequences.txt"
	shell:"""
		orthomclBlastParser {input.blastfile} compliantFasta/ >> {output}
	"""
	

#Since LoadBlast doesn't generate a file, I create a control (empty) one. 
rule LoadBlast:
	input:
		"similarSequences.txt"
	output:
		"controlLoadBlast"
	shell:"""
		orthomclLoadBlast orthomcl.config {input};
		touch controlLoadBlast
	"""

rule Pairs:
	input:
		"controlLoadBlast"
	output:
		"ortho10.log"
	shell:"""
		orthomclPairs orthomcl.config ortho10.log cleanup=no
	"""

rule DumpPairsFiles:
	input:
		"ortho10.log"
	output:
		"mclInput"
	shell:"""
		orthomclDumpPairsFiles orthomcl.config
	"""

rule mcl:
	input:
		"mclInput"
	output:
		"mclOutput"
	shell:"""
		mcl {input} --abc -I 1.5 -o {output}
	"""

rule mclToGroups:
	input:
		"mclOutput"
	output:
		"groups.txt"
	shell:"""
		orthomclMclToGroups micrurus 0000000 < {input} > {output}	
	"""
