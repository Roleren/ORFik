ORFik: R package for discovery of novel genes.
==============================================================================

This package is still under development, although this version is stable and can be used.

#### About


ORFik is an R package containing various functions for analysis of Ribo-Seq, RNA-Seq and Cage-Seq data. ORFik currently supports:

- Finding Open Reading Frames (very fast) in the genome of interest or on the set of transcripts.
- Metaplots for Ribo-Seq allowing to spot the shift.
- Automatically detecting shifts for Ribo-Seq data.
- Various measurements of gene identity e.g. FLOSS, coverage, ORFscore, entropy.
- Utility functions to extend GenomicRanges, e.g. the functions tile1 and groupGRangesBy.



#### Installation
Package is available from bioconductor (3.8, R version >= 3.5.0)
```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ORFik")
```

Development version on bioconductor (3.9, R version >= 3.6.0)
```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ORFik", version = "devel")
```  

Package is also available here on github
```r
library(devtools)
install_github("JokingHero/ORFik")
```  

#### More information

After installation run:
```r
library(ORFik)

# Tissue specific 5' utrs using cage-data
?reassignTSSbyCage

# Detecting open reading frames
?findORFs

# get a feature-set from predicted orfs
?computeFeatures

# read vignette
browseVignettes("ORFik")
```  

#### Feedback

Please feel free to provide feedback or desired functionality. My contact address is kornellabun@gmail.com.
