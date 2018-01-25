ORFik: R package for discovery of novel genes with Ribo-Seq and RNA-Seq data.
==============================================================================

This package is still under development, although this version is stable and can be used already.

#### About


ORFik is an R package containing various functions for analysis of Ribo-Seq and RNA-Seq data. Currently it supports:

- Finding Open Reading Frames (very fast) in the genome of interest or on the set of transcripts.
- Metaplots for Ribo-Seq allowing to spot the shift.
- Shifting functions for the Ribo-Seq data.
- Various measurements of gene identity eg. FLOSS, coverage, ORFscore, entropy.



#### Installation
Package is currently only available on github
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
?find_in_frame_ORFs

# get a feature-set from predicted orfs
?allFeatures

# to be added:
# By machine learning, rank the predicted
# orfs, from the the feature set.

# read vignette
browseVignettes("ORFik")
```  

#### Feedback

Please feel free to provide feedback or desired functionality. My contact address is Kornel.Labun at uib.no.
