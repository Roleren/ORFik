ORFik: R package for discovery of novel genes.
==============================================================================

This package is still under development, although this version is stable and can be used.

#### About


ORFik is a R package containing various functions for analysis of Ribo-Seq, RNA-Seq, Cage-Seq and TCP-seq data related to transcriptomics. ORFik currently supports:

- Finding Open Reading Frames (very fast) in the genome of interest or on the set of transcripts.
- Metaplots for Ribo-Seq allowing to spot the shift.
- Automatically detecting shifts for Ribo-Seq data.
- Various measurements of gene identity e.g. FLOSS, coverage, ORFscore, entropy.
- CAGE reassignment of leaders in txdb/GenomicRanges objects.
- Easy creation of coverage plots of sequence data over transcript regions.
- Utility functions to extend and speed up GenomicRanges, e.g. the functions tile1, groupGRangesBy and pmapFromTranscriptF.



#### Installation
Package is available from bioconductor (3.9, R version >= 3.6.0)
```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ORFik")
```

Development version on bioconductor (3.10, R version >= 3.6.0)
```r
if (!requireNamespace("BiocManager", quietly=TRUE))
    install.packages("BiocManager")
BiocManager::install("ORFik", version = "devel")
```  

Package is also available here on github
```r
library(devtools)
install_github("Roleren/ORFik")
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

Please feel free to provide feedback or desired functionality by creating a new issue on our github page.
