#####################################################################################
# groups2fasta.py                                                                   #
# Generate the fasta files containing the family genes from the ortholog groups.    #
# input: groups.txt from orthomcl output AND fasta CDS files for each sample.       #
# output: fasta files with the ortholog RNA sequences.                              #
#                                                                                   #
## Considerations:                                                                  #
# - The fasta files should be renamed with the same 4 letters codes used            #
# in the previous orhomcl run.                                                      #
#                                                                                   #
#                                                                                   #
## Example run:                                                                     #
# python groups2fasta.py -j 12  --config fastadir=../data/ -np                      #
#####################################################################################

import subprocess
from os.path import join
import os,re,sys
from Bio import SeqIO
from os.path import basename


# Globals ---------------------------------------------------------------------

groups="groups.txt"

# Full path to a folder that holds all of your FASTA files.
FASTA_DIR = config["fastadir"]

# A Snakemake regular expression matching the FASTA files.
SAMPLES, = glob_wildcards(join(FASTA_DIR, '{sample,[^/]+}.fasta'))

# Patterns using the 'sample' wildcard.
PATTERN = '{sample}.fasta'

# Functions -------------------------------------------------------------------

# Rules -----------------------------------------------------------------------

rule all:
	input:
		"CDSgroups/"

rule adjustFASTA:
	input:
		files=join(FASTA_DIR, PATTERN),
	output:
		compliantFiles='CDScompliantFasta/{sample}.fasta'
	shell:"""
		name=`echo -n {input.files} | tail -c10 | head -c4`;
		orthomclAdjustFasta $name {input.files} 2;
		mv $name.fasta CDScompliantFasta/;
	"""

#Keep only the clusters with reads from all the distinct input samples.
rule filterGroups:
	input: 
		compliantFiles=expand('CDScompliantFasta/{sample}.fasta', sample=SAMPLES),
		groupfile=groups
	output:
		'groups_OK.txt'
	run:
		lineok=""
		for line in open(groups):
			recordlist=[]
			for index,record in enumerate(line.split("|")[:-1]):
				recordlist.append(record[-4:])
			if len(list(set(recordlist)))==len(input.compliantFiles):
				lineok+=line
		with open(output[0],'w') as f:
			f.write(lineok)

rule extractGroups:
	input:
		'groups_OK.txt'
	output:
		'CDSgroups/'
	run:
		aux=[]
		nameFile=''
		with open(input[0]) as f:
		    for line in f:
		        file_sequences = []
		        aux=line.split(': ')
		        nameFile=aux[0]
		        path="CDScompliantFasta"
		        aux=aux[1].split(' ')
		        for id in aux:
		            gen_rec=id.split("|")
		            output_res = open("CDSgroups/"+nameFile+'.fa', 'w')
		            for record in SeqIO.parse(open(path+'/'+gen_rec[0]+'.fasta', "rU"), "fasta") :
		                #print "A: "+record.id
		                #print "B: "+gen_rec[0]+"|"+gen_rec[1].rstrip()
		                if record.id == gen_rec[0]+"|"+gen_rec[1].rstrip():
		                    #print "fasta id: "+record.id
		                    #print "Lista id: "+gen_rec[0]+"|"+gen_rec[1].rstrip()+"\n"
		                    #record.id=gen_rec[0]+"|"
		                    file_sequences.append(record)
		                    break
		        SeqIO.write(file_sequences, output_res, "fasta")
