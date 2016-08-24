###########################################################################################################
### From a batch of fasta files of family genes (orthomcl groups output) it obtains a phylogenetic tree.  #
# Steps:                                                                                                  #
# 1. Remove duplicate records. One record from each sample is kept (longest).                             #
# 2. Convert clusters to protein and do alignment with muscle.                                            #
# 3. Arrange alignments by ID.                                                                            #
# 4. Concatenate and create partitions file with length coordinates.                                      #
# 5. Convert fasta concatenated file to phy.                                                              #
# 6. Build phylogenetic tree: raxml.                                                                      #
                                                                                                          #
### Considerations:                                                                                       #
# - All the fasta files should contain records (at least one) from ALL the samples in our study.          #
# - The IDs from fasta files should have only 4 letters (orthomcl requeriments)(>XXXX|).                  #
# - For translation to protein, transeq or biopython are available (default biopython).                   #  
# - For remove duplicates, two options: remove the shorter ones or remove them randomly (default random). #
# - Default raxml bootstrap: 100.                                                                         # 
#                                                                                                         #
###Prerequisites                                                                                          #
# - transeq from EMBOSS, in case the user prefer to use it (instead biopython) to perform                 #
# the protein translation.                                                                                #
# - Biotpython and bioawk.                                                                                #
                                                                                                          #
### Usage (slurm example):                                                                                #                                                                            #
# $ srun --partition=compute --time 7-0 --mem=10G --cpus-per-task=24 --ntasks=1 --pty bash                #
# $ module load python/3.5.0                                                                              #
# (dry run) $ snakemake --snakefile snakephylo -j 24 --config fastadir=data/ -np                          #
# $ snakemake --snakefile snakephylo -j 24 --config fastadir=data/                                        #
###########################################################################################################

import subprocess,sys
from os.path import join
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.SeqRecord import SeqRecord

# Globals ---------------------------------------------------------------------

# Full path to a folder that holds all of your FASTA files.
#FASTA_DIR = 'in/'
FASTA_DIR = config["fastadir"]

# A Snakemake regular expression matching the FASTA files.
SAMPLES, = glob_wildcards(join(FASTA_DIR, '{sample,[^/]+}.fa'))

# Patterns using the 'sample' wildcard.
PATTERN = '{sample}.fa'

# Functions -------------------------------------------------------------------
#
def sortAlign(sample):
    res=""
    handle = open(str(sample), "rU")
    l = SeqIO.parse(handle, "fasta")
    sortedList = [f for f in sorted(l, key=lambda x : x.id)]
    for s in sortedList:
        res+=">"+s.description+"\n"+str(s.seq)+"\n"
    return res

#Check repetated records from a sample file and delete them randomly
def deleteRepeated(filename):
    lista=[]
    res=""
    with open(str(filename), 'r') as f:
        for line in f:
            if line[0]==">":
                id=line[1:5]
                if id not in lista:
                    lista.append(id)
                    copy=1
                else: copy=0
            if copy: res+=line
    return res

#Check repetated records from a sample file and delete the shorter ones
def keepLongerRead(filename):
    fasta_sequences = SeqIO.parse(open(str(filename), "rU"),'fasta')
    recdict={}
    for seq_record in fasta_sequences:
        i= seq_record.id.split("|")[0]
        if i in recdict:
            if len(seq_record)>len(recdict[i]):
                recdict[i]=[seq_record,seq_record.description]
        else:
            recdict[i]=[seq_record,seq_record.description]
    records=[]
    for value in recdict:
        records.append(SeqRecord(Seq(str(recdict[value][0].seq), generic_dna), id=value+"|"+recdict[value][1][5:], description=''))
    return records

def make_protein_record(nuc_record):
    """Returns a new SeqRecord with the translated sequence (default table)."""
    return SeqRecord(seq = nuc_record.seq.translate(cds=False),id = nuc_record.id, description = "")

# Rules -----------------------------------------------------------------------
#
rule all:
    input:
        'RAxML_bestTree.outRaxml'

rule removeDuplicates:
    input:
        join(FASTA_DIR, PATTERN)
    output:
        'nodup/{sample}.fasta'
    run:
        with open(output[0],'w') as f:
            f.write(deleteRepeated(input))
        #output_handle = open(output[0], "w")
        #SeqIO.write(keepLongerRead(input), output_handle, "fasta")
           
rule toProtein:
    input:
        'nodup/{sample}.fasta'
    output:
        'prot/{sample}.fasta'
    #2 options, biopython or transeq (EMBOSS)
    run:
        proteins = (make_protein_record(nuc_rec) for nuc_rec in SeqIO.parse(str(input), "fasta"))
        SeqIO.write(proteins, output[0], "fasta")
    #shell:"""
    #    transeq --sequence {input} -outseq {output} -sformat pearson
    #"""

rule align_protein:
    input:
        'prot/{sample}.fasta'
    output:
        'muscled/{sample}.fasta'
    shell:"""
        muscle -in {input} -out {output}
    """

rule sort_alignment:
    input:
        'muscled/{sample}.fasta'
    output:
        'muscled/sorted/{sample}.fasta'
    run:
        with open(output[0],'w') as f:
            f.write(sortAlign(input))            

rule concatenate:
    input: 
        expand('muscled/sorted/{sample}.fasta', sample=SAMPLES)
    output:
        concat='concat.fasta',
        parti='partitions.csv'
    run:
        numlines = sum(1 for line in open(input[0]))
        mixed = [""] * int(numlines)
        indexRecord=1
        res=""
        for indexFile,i in enumerate(input):
            lenRecord = subprocess.check_output("bioawk -c fastx '{ print $name, length($seq) }' < "+i+" | awk '{print $2}' | uniq", shell=True)
            res+="WAG, "+ i.split("/")[-1][:-6] +" = "+str(indexRecord)+"-"+str(int(lenRecord[:-1])+indexRecord-1)+"\n"
            indexRecord=int(lenRecord[:-1])+indexRecord
            for index,line in enumerate(open(i,'r')):
                if indexFile==0 or (indexFile!=0 and index%2):
                    if line[0]==">": mixed[index]+=line[:5]
                    else: mixed[index]+=line[:-1]
        with open(output.concat,'w') as f:
            for item in mixed:
                print(item,sep="\n",end="\n",file=f)
        with open(output.parti, 'w') as out:
            out.write(res)

rule toPHY:
    input:
        'concat.fasta'
    output:
        'concat.phy'
    shell:"""
        prank -convert -d={input} -f=phylipi -o=concat -keep
    """

rule raxml:
    input:
        partitions = 'partitions.csv',
        phy = 'concat.phy' 
    output:
        'RAxML_bestTree.outRaxml'
    shell:"""
        module load raxml/8.2.4
        raxml -T 24 -m PROTGAMMAWAG -s {input.phy} -n outRaxml -p 12345 -q {input.partitions} -x 12345 -f a -# 100
    """    
#without bootstrap
#raxml -T 24 -m PROTGAMMAWAG -s {input.phy} -n outRaxml -p 12345 -q {input.partitions}
