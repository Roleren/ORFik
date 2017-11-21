ORFik: R package for finding Open reading frames and creating ORF features
==============================================================================

This package is still under development, although this version is stable and can be used already.

#### About
Our R package’s scans for Open Reading Frames, giving GenomicRanges compatitable output.
Functions for the most popular features used are implemented, like: coverage, entropy etc.
Meta plots for easy visualization, are also provided.

#### Installation
Package is only available on github
```r
library(devtools)
install_github("JokingHero/ORFik")
```  

#### More information

After installation run:
```r
library(ORFik)

# main function
?find_in_frame_ORFs

# read vignette
# Not implemented yet
#browseVignettes("ORFik")
```  

#### Feedback

Please feel free to provide feedback or desired functionality. My contact address is Kornel.Labun at uib.no.
