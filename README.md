# Motif Mark

## Objective

The motifs, or small regions of the DNA that are shared upon different different proteins, can have an impact on the important biological functions of a gene. Multiple motifs can work together and regulate gene expression. 

The objective of this motif-mark is to determine and visualize motif-binding site on DNA sequences.

##

## Setting up environment
The python script uses Pycairo to create a figure with motifs. Python 3 and Pycairo (version) should be installed in the environment.

## How to use
In the working directory, there should be 3 files:

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

1. FASTA file: The fasta file contains a header and DNA sequences of the reads.
2. Motif text file: The file contains a motif on every line. This script can accept ambiguous nucleotide motifs.

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

1. oneline_fasta.txt : This file contains the FASTA file, concatenating the multiple lines of sequence into one sequence line. Therefore, this file alternates between header and sequence line for every new line.
2. 'fastafilename'.png : The png will have the same file name as the FASTA file. The png will have visualizations of motifs for every reads.

## Things to work on:
Motif colors are hardcoded to only have 4 colors: purple, pink, yellow, and blue. If there are more motifs in the motif file, not all motifs will be visualized in the figure. To fix this, I will create a random number generator to create random numbers. The exact same 3 numbers will not be used to ensure no same color will be accidently created.

The python script does not take all ambiguous nucleotides. Only Y. Ensure it can take all the ambivuous nucleotides, including the above.
