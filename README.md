# Fu2023

Repository for files related to the academic paper "Identification of genes with oscillatory expression in glioblastoma – The paradigm of SOX2" by Fu et. al 2023. We have chosen to split the repository into 3 separate sub-directories so each can be reproduced individually. You need only use the codes in the sub-directory relevant to the specific pipeline you would like to perform. Detailed instructions on how to install and use these pipelines can be found in READMEs within the relevant sub-directory. 

## Network_Screen 

Inference of transcription factor and microRNA regulatory interactions using ChIP-seq and miRNA target data. This was used in this study to predict oscillatory transcription factors by virtue of their involvement in dynamic(al) gene network motifs. This was used for the analysis in Figure 2c. 

This network of transcriptional interactions with miRNA regulation was first conceived and constructed by Tom Minchington as part of his PhD in the Nancy Papalopulu lab. This work was published and can be found here (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7395740/).

The scripts in this repository can be used to build an updated version of this network, using the latest databases available as of May 2023. This new network was generated as part of the experimental workflow detailed in R. Fu et al 2023 by Oliver Cottrell.

## OscoNet

Algorithm to infer oscillatory genes from scRNA-seq data. This was performed by Richard Fu to generate the OscoNet data in Figure 2 and beyond.

This pipeline was first conceived and published in Cutillo et al., 2020 (https://doi.org/10.1186/s12859-020-03561-y). Here we provide a simplified, more user-friendly distribution of the original method, which was modified by Oliver Cottrell. No changes have been made to the core functionality of the code.

## Periodicity_Analysis

Matlab pipeline to determine periodicity from SOX2-mKate2 GBM1 cells. This code was writted by Andrew Rowntree to construct Figure 6 of Fu et al., 2023.

