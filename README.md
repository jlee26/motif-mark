# Motif Mark

## Objective

The motifs, or small regions of the DNA that are shared between different different proteins, can have an impact on biological functions of a gene. Multiple motif binding sites can work together and regulate gene expression. These sequences can be as short as 3 nucleotides and as along as over 10 nucleotides. It is possible some motif binding sites are regulators of alternative splicing, a biological process where introns are spliced and exons are joined together for protein synthesis, where certain motif binding sites can be an activator or repressor for exon inclusion.

The objective of this motif-mark script is to determine and visualize motif-binding sites on DNA sequences. The script will return one figure for a single FASTA file. In the figure, motif binding sites will be colored in boxes to scale with the DNA sequence. This figure is to help understand the location of the motif binding sites and support biological conclusions.

## Setting up environment
The motif-mark script uses Pycairo to create a figure with motif binding sites. Python (3.8.12) and Pycairo should be installed to the environment.

## How to use
In the working directory, there should be 3 files (more described below):

1. motif-mark-oop.py
2. FASTA file
3. Motif text file

The following options are required to use the python script:
- -f --file : FASTA file
- -m --motif : Motif file

Type the following in the command line:
```
./motif-mark-oop.py -f fastafile.fasta -m motif.txt
```

Example:
```
./motif-mark-oop.py -f Figure_1.fasta -m Fig_1_motifs.txt
```

## Input
The python script will require two files:

1. FASTA file: The fasta file contains a header and DNA sequences of the reads. The DNA sequences are case sensitive. Upper cases indicate exons and lower cases indicate introns.
2. Motif text file: The file contains a motif on every line. This script can accept ambiguous nucleotide motifs. This file is NOT case sensitive.

| Description  |     Symbol    | Bases represented |
| ------------ | ------------- | ----------------- |
|   Adenine    |       A       |          A        |
|   Cytosine   |       C       |          C        |
|   Guanine    |       G       |          G        |
|   Thymine    |       T       |          T        |
|    Uracil    |       U       |        U/T        |
|     Weak     |       W       |        A/T        |
|    Strong    |       S       |        C/G        |
|     Amino    |       M       |        A/C        |
|    Ketone    |       K       |        G/T        |
|    Purine    |       R       |        A/G        |
|  Pyrimidine  |       Y       |        C/T        |
| Any one base |       N       |     A/C/G/T       |

## Output
The python script will return 2 output files.

1. oneline_fasta.txt : This file is the FASTA file concatenating the multiple lines of sequence into one sequence line.
2. 'fastafilename'.png : The png will have the same file name as the FASTA file. The png will have visualizations of motifs for every reads. It is possible the motif binding site colors are similar and hard to distinguish from other motifs. It is highly recommended to rerun the script as many times to result in contrasting colors.

## Future Work
The figure will automatically create different colors for each motif binding sites. The script currently does not support similar colors to not be used. If figure has too many similar colored motifs, rerun the script to result in more contrasting colors.
