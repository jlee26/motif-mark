# Motif Mark

## Objective

The motifs, or small regions of the DNA that are shared upon different different proteins, can have an impact on the important biological functions of a gene. Multiple motifs can work together and regulate gene expression. 

The objective of this motif-mark is to determine and visualize motif-binding site on DNA sequences.

##

## Setting up environment
The python script uses Pycairo to create a figure with motifs. Python 3 and Pycairo (version) should be installed in the environment.

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
