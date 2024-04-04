ORFik: R package for discovery of novel genes.
==============================================================================

[![bioc](http://www.bioconductor.org/shields/years-in-bioc/ORFik.svg)](http://bioconductor.org/packages/devel/bioc/html/ORFik.html) [![bioc](http://www.bioconductor.org/shields/downloads/devel/ORFik.svg)](https://bioconductor.org/packages/stats/bioc/ORFik/)  <a href="https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-021-04254-w"><img border="0" src="https://img.shields.io/badge/published-BMC-blue.svg" title="Published in BMC bioinformatics"></a>

![](inst/images/ORFik_map.png)

This package is published, but still under heavy development to include more features.

#### About


ORFik is a R package containing various functions for analysis of Ribo-Seq, RNA-Seq, CAGE and TCP-seq data related to transcriptomics. ORFik currently supports:

1. Finding Open Reading Frames (very fast) in the genome of interest or on the 
set of transcripts/sequences.  
2. Hundreds of functions helping your analysis of either: sequence data, RNA-seq data, CAGE data, Ribo-seq data, TCP-seq data or RCP-seq data.
3. Automatic estimations of RiboSeq footprint shift.  
4. Utilities for metaplots of RiboSeq coverage over gene START and STOP codons 
allowing to spot the shift.  
5. Shifting functions for the RiboSeq data.  
6. Annotation / re-annotation of 5' UTR Transcription Start Sites using CAGE data.  
7. Various measurements of gene identity, more than 30 functions. e.g. FLOSS, coverage, ORFscore, 
entropy that are recreated based on scientific publications.  
8. Utility functions to extend GenomicRanges for faster grouping, splitting, filtering etc. Included c++ function for speed.
9. Extensive implemented syntax for coverage and metacoverage of NGS data, including smart grouping functions for easier prototyping.
10. Automatic download of genome annotation from any species supported by ensembl.
11. Automatic download and metadata extraction of NGS files from SRA, ERA, DRA and GEO. 
12. Full NGS alignment pipeline: Trimming data using fastp and alignment using STAR (with optional contaminant removals)
13. Simplifying working with massive amounts of datasets using the ORFik experiment class. 


#### Installation
Package is available from bioconductor (Stable branch, R version >= 4.0.0)
```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ORFik")
```

Development version on bioconductor (Devel branch, R version >= 4.0.0)
```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ORFik", version = "devel")
```  

Package is also available here on github (Experimental branch, R version >= 4.1.0)
```r
if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
remotes::install_github("Roleren/ORFik")
```  

#### More information

After installation run:
```r
library(ORFik)

# NGS metadata extraction
?download.SRA.metadata

# NGS data download
?download.SRA

# Annotation download
?getGenomeAndAnnotation

# Data management
?create.experiment

# NGS Library Quality control
?ORFikQC

# Tissue specific 5' utrs using cage-data
?reassignTSSbyCage

# Detecting open reading frames
?findORFs

# get a feature-set from predicted orfs
?computeFeatures

# read vignette (tutorials)
browseVignettes("ORFik")
```  
Please read Bioconductor vignettes for detailed tutorials and examples.

#### Feedback

Please feel free to provide feedback or desired functionality by creating a new issue on our github page.
