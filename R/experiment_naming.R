#' Get library name variants
#'
#' Used to standardize nomeclature for experiments.\cr
#' Example: RFP is main naming, but a variant is ribo-seq
#' ribo-seq will then be renamed to RFP
#' @family experiment_naming
#' @return a data.table with 2 columns, the main name, and all name variants
#' of the main name in second column as a list.
libNames <- function() {
  mainName <- c("RNA", "RFP", "CAGE", "LSU",
                "SSU", "ATAC", "tRNA", "SHAPE",
                "ChIP", "PRPF")
  allNames <-
    list(c("rna-seq", "Rna-seq", "RNA-seq", "RNA-Seq"),
         c("RFP", "RPF", "ribo-seq", "Ribo-seq", "ribo-Seq"),
         c("CAGE", "cage"),
         c("80S","LSU"),
         c("40S","SSU"),
         c("ATAC"),
         c("PRPF"),
         c("SHAPE"),
         c("ChIP"),
         c("tRNA")
    )
  dt <- data.table(mainName, allNames)
  return(dt)
}

#' Get stage name variants
#'
#' Used to standardize nomeclature for experiments.\cr
#' Example: 64Cell stage is same as 2 hours post fertilization,
#' so all 2hpf will be converted to 64Cell etc.
#' @family experiment_naming
#' @return a data.table with 2 columns, the main name, and all name variants
#' of the main name in second column as a list.
stageNames <- function() {
  mainName <- c("unfertilized", "fertilized",
                "2to4Cell", "4Cell", "8Cell", "64Cell", "256Cell", "512Cell",
                "1KCell", "High", "Oblong", "Sphere", "Dome", "Shield", "Bud",
                "Somite", "24hpf", "prim6", "prim10", "prim12",  "prim20",
                "2dpf", "3dpf", "4dpf", "5dpf", "6dpf", "10dpf", "21dpf",
                "24dpf")
  allNames <-
    list(c("unfertilized", "Unfertilized"),
         c("_fertilized", "_Fertilized"),
         c("2to4Cell", "2to4cell", "2to4_cell", "2-4cell", "2-4Cell", "2-4_cell"),
         c("4cell", "4Cell", "4_cell"),
         c("8cell", "8Cell", "8_cell"),
         c("64cell", "64Cell", "64_cell", "64_Cell", "_2h", "_02h", "2hpf"),
         c("256cell", "256Cell", "256_cell", "256_Cell"),
         c("512cell", "512Cell", "512_cell"),
         c("1Kcell", "1KCell", "1K_cell", "_3h", "_03h", "3hpf"),
         c("High", "high"),
         c("Oblong", "oblong"),
         c("Sphere", "sphere", "_4h", "_04h","4hpf"),
         c("Dome", "dome"),
         c("Shield", "shield", "_6h", "_06h", "6hpf"),
         c("Bud", "bud", "_10h", "10hpf"),
         c("Somite", "somite", "_12h", "12hpf"),
         c("24hpf", "_24h", "1dfp"),
         c("prim6", "prim_6", "25hpf", "_25h"),
         c("prim10", "prim_10", "27hpf", "_27h"),
         c("prim12", "prim_12", "28hpf", "_28h"),
         c("prim20", "prim_20", "33hpf", "_33h"),
         c("2dpf", "_48h", "_48hpf"),
         "3dpf", "4dpf", "5dpf", "6dpf", "10dpf", "21dpf", "24dpf"
         )
  dt <- data.table(mainName, allNames)
  return(dt)
}

#' Get tissue name variants
#'
#' Used to standardize nomeclature for experiments.\cr
#' Example: testis is main naming, but a variant is testicles.
#' testicles will then be renamed to testis.
#' @family experiment_naming
#' @return a data.table with 2 columns, the main name, and all name variants
#' of the main name in second column as a list.
tissueNames <- function() {
  mainName <- c("adipose", "amygdala","brain", "bladder", "blood", "bone","breast",
    "colon", "cortex", "eye", "fibroblast", "frontal_lobe", "heart",
    "kidney", "liver", "lung", "melanocyte", "mesenchymal", "monocyte", "muscle",
    "myeloid","ovary", "prostate",
    "rectum", "retina","testis","urunary", "vagina", "skin", "tongue")
  allNames <-
    list(c("adipose", "Adipose", "Adipocyte", "adipocyte"),
         c("amygdala", "Amygdala"),
         c("brain", "Brain"),
         c("bladder", "Bladder"),
         c("blood", "Blood"),
         c("bone", "Bone", "osteosarcoma", "osteogenic"),
         c("breast", "Breast"),
         c("colon", "Colon"),
         c("cortex", "Cortex"),
         c("eye", "Eye"),
         c("fibroblast", "Fibroblast"),
         c("frontal_lobe", "frontal lobe", "Frontal_lobe"),
         c("heart", "Heart"),
         c("kidney", "Kidney", "Renal", "renal"),
         c("liver", "Liver"),
         c("lung", "Lung"),
         c("Melanocyte", "melanocyte"),
         c("mesenchym", "Mesenchym"),
         c("Monocyte", "monocyte"),
         c("muscle", "Muscle", "myotubes", "Myotubes"),
         c("myeloid", "Myeloid"),
         c("ovary", "Ovary"),
         c("prostate", "Prostate"),
         c("rectum", "Rectum", "anus", "Anus"),
         c("Retina", "retina"),
         c("testis", "Testis", "testicles", "Testicles"),
         c("urunary", "Urunary"),
         c("vagina", "Vagina"),
         c("skin", "Skin"),
         c("tongue", "Tongue")
    )

  dt <- data.table(mainName, allNames)
  return(dt)
}
#' Get cell-line name variants
#'
#' Used to standardize nomeclature for experiments.\cr
#' Example: THP-1 is main naming, but a variant is THP1
#' THP1 will then be renamed to THP-1
#' @family experiment_naming
#' @return a data.table with 2 columns, the main name, and all name variants
#' of the main name in second column as a list.
cellLineNames <- function() {
  mainName <- c("HEK293", "HeLa", "THP1", "PC3")
  allNames <-
    list(c("HEK293", "Hek293"),
         c("HeLa", "HELA", "Hela", "hela"),
         c("THP-1", "THP1"),
         c("PC3", "PC-3")
    )
  dt <- data.table(mainName, allNames)
  return(dt)
}

#' Get replicate name variants
#'
#' Used to standardize nomeclature for experiments.\cr
#' Example: 1 is main naming, but a variant is rep1
#' rep1 will then be renamed to 1
#' @family experiment_naming
#' @return a data.table with 2 columns, the main name, and all name variants
#' of the main name in second column as a list.
repNames <- function() {
  mainName <- c("1", "2", "3", "4", "5", "6")
  allNames <-
    list(c("rep1", "Rep1", "run1", "_r1_"),
         c("rep2", "Rep2", "run2", "_r2_"),
         c("rep3", "Rep3", "run3", "_r3_"),
         c("rep4", "Rep4", "run4", "_r4_"),
         c("rep5", "Rep5", "run5", "_r5_"),
         c("rep6", "Rep6", "run6", "_r6_")
    )
  dt <- data.table(mainName, allNames)
  return(dt)
}

#' Get condition name variants
#'
#' Used to standardize nomeclature for experiments.\cr
#' Example: WT is main naming, but a variant is control
#' control will then be renamed to WT
#' @family experiment_naming
#' @return a data.table with 2 columns, the main name, and all name variants
#' of the main name in second column as a list.
conditionNames <- function() {
  mainName <- c("WT", "MZ", "4Ei", "Silvesterol",
                "Mutant", "cas9", "NMDA", "DHPG")
  allNames <-
    list(c("WT", "wt", "wild_type", "Wild_type",
           "control", "Control", "Basal"),
         c("MZ", "dicer"),
         c("4Ei", "4ei"),
         c("silvesterol", "Silvesterol"),
         c("mutant", "Mutant"),
         c("Cas9", "cas9"),
         c("NMDA"),
         c("DHPG")
    )
  dt <- data.table(mainName, allNames)
  return(dt)
}

#' Get main name from variant name
#'
#' Used to standardize nomeclature for experiments.\cr
#' Example: RFP is main naming, but a variant is ribo-seq
#' ribo-seq will then be renamed to RFP
#' @param names a character vector of names that must exist in dt$allNames
#' @param dt a data.table with 2 columns (mainName, allNames)
#' @family experiment_naming
#' @return a data.table with 2 columns, the main name, and all name variants
#' of the main name in second column as a list.
mainNames <- function(names, dt) {
  g <- groupings(dt$allNames)
  dt.long <- data.table(mainName = c(dt$mainName[g], c("", "")),
                        allNames = c(unlist(dt$allNames), c("", NA)))
  return(dt.long[chmatch(names, dt.long$allNames),]$mainName)
}
