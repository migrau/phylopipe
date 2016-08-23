# Phylogeny pipeline

## Description

Scripts used to perform a phylogeny.

## Contents

* _src/snakeorthomcl.py_. Run orthomcl using a batch of protein fasta files.
* _src/snakephylo.py_. Obtain raxml tree from a batch of unaligned fasta files, using codon mode.
* _data_. Required fasta files to complete an example run.

### snakeorthomcl

_OrthoMCL Snakemake script_

It runs all the steps of orthoMCL automatically, parallelizing blast (based on http://orthomcl.org/common/downloads/software/v2.0/UserGuide.txt).
- Input files: CDS proteins (obtanied from transdecoder for example).
- Output file: groups.txt, containing all the ortholog clusters.

_Considerations_

- An empty mysql database is required.
- FASTA files should be renamed using a code of only four letters: XXXX.fasta
- orthomcl.config file is also required (attached).
- Blast step is parallelized, by default 50 threads.

_Possible error_

1. For rule adjustFASTA, the user has to indicate the position of the unique ID. By default, it is 2, change if necessary.

_Usage_

```{bash}
$ snakemake --snakefile snakeorthomcl.py (-np) -j 24
```

### snakephylo

From a batch of fasta files of family genes (orthomcl groups output) it obtains a phylogenetic tree.  

_Steps_

1. Remove duplicate records. One record from each sample is kept (longest or random).
2. Translate clusters to protein and do alignment with muscle.
3. Arrange alignments by ID.
4. Concatenate and create partitions file with length coordinates.
5. Convert fasta concatenated file to phy.
6. Build phylogenetic tree: raxml.

_Considerations_ 

- All the fasta files should contain records (at least one) from ALL the samples in our study. 
- The IDs from fasta files should have only 4 letters (orthomcl requeriments)(>XXXX|). 
- For translation to protein, transeq or biopython are available (default biopython).
- For remove duplicates, two options: remove the shorter ones or remove them randomly (default random).
- Default raxml bootstrap: 100. 
 
_Prerequisites_

- transeq from EMBOSS, in case the user prefer to use it (instead biopython) to perform the protein translation.  
- Biotpython and bioawk are required in any case. 

_Usage (slurm example)_

```{bash}
$ srun --partition=compute --time 7-0 --mem=10G --cpus-per-task=24 --ntasks=1 --pty bash 
$ module load python/3.5.0 
(dry run) $ snakemake --snakefile snakephylo.py -j 24 --config fastadir=data/ -np 
$ snakemake --snakefile snakephylo.py -j 24 --config fastadir=data/ 
```

## Steps
