## RNA-seq Differential Expression and Pathway Analysis  
This project involved a complete RNA-seq analysis workflow to identify and interpret transcriptional changes between control and treated samples in human T cells.

## Overview  
The goal of this analysis was to reproduce and extend the original findings from an RNA-seq study examining the effect of IL-27 on T cell activation. The pipeline includes preprocessing, quality control, alignment, quantification, normalization, differential gene expression, and pathway enrichment analysis.

## Objectives  
- Build a complete RNA-seq analysis pipeline using R and shell scripting.  
- Identify differentially expressed genes (DEGs) between control and treatment groups.  
- Interpret enriched biological pathways associated with IL-27 stimulation.  
- Compare findings with the original study results and highlight novel insights.  

## Tools & Technologies  
- Tools: FastQC, STAR, SAMtools, featureCounts, DESeq2, clusterProfiler.  
- Data Sources: Public RNA-seq datasets (ENA), gene annotation from Ensembl.  
- Environment: Linux, command-line interface.
  
  ## Pipeline Components
1. **Sample Download**  
   Automated download of raw FASTQ files from ENA.

2. **Quality Control**  
   Initial assessment of read quality using FastQC.

3. **Genome Indexing & Alignment**  
   - STAR used to index the genome and align reads  
   - SAMtools used to sort and convert alignment files

4. **Quantification**  
   - featureCounts for generating raw count matrices  
   - Scripts to organize and merge count data for DESeq2

5. **Normalization & DEG Analysis**  
   - DESeq2 used for normalization and statistical testing  
   - Generation of MA plots, PCA plots, and volcano plots

6. **Pathway Analysis**  
   - clusterProfiler used for GO and KEGG enrichment  
   - Focused interpretation of immune regulation pathways

> **Data Source**: This project reanalyzes RNA-seq data from the study _"Leveraging transcription factors to speed cellobiose fermentation by Saccharomyces cerevisiae"_  
> [View the original publication](https://experts.illinois.edu/en/publications/leveraging-transcription-factors-to-speed-cellobiose-fermentation)
