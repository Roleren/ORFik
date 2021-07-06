#' Get library name variants
#'
#' Used to standardize nomeclature for experiments.\cr
#' Example: RFP is main naming, but a variant is ribo-seq
#' ribo-seq will then be renamed to RFP
#' @family experiment_naming
#' @return a data.table with 2 columns, the main name, and all name variants
#' of the main name in second column as a list.
libNames <- function() {
  mainName <- c("RNA", "RFP", "QTI", "CAGE", "LSU",
                "SSU", "ATAC", "PRPF", "PAS-Seq", "PAL-Seq", "RIP","SHAPE",
                "ChIP", "CLIP", "tRNA", "miRNA", "GRO")
  allNames <-
    list(c("rna-seq", "Rna-seq", "RNA-seq", "RNA-Seq", "RNASeq", "RNAseq", "RNA seq", "RNASEQ",
           "rnaseq", "Input", "input", "total RNA", "Total RNA", "total_RNA", "totalRNA","TotalRNA_", "_Total_",
           "mRNA$", " mRNA ","_mrna_", "_RNA_", "_rna_", "^rna_", "^RNA_", "^mrna_", "^RNA ", "^mRNA ", "^mRNA_",
           "\\.rna$", "\\.mrna", "_RNA$", "_rna$","PolyA"),
         c("RFP", "RPF", "ribo-seq", "ribo-Seq", "Ribo-seq", "Ribo-Seq", "Ribo_seq",
           "riboseq", "Riboseq", "RiboSeq", "Ribosome", "ribosome", "Profiling","profiling", "Ribo\\.prof",
           "Footprint", "ribo_profile", "footprint", "^Ribo_", "^ribo_", "^Ribo ","^ribo_", "_ribo_", "_Ribo_",
           "_ribo$", "_Ribo$", "\\.ribo",
           "^RP ", " RP ", "_RP_", "^FP_", "^fp_", "_fp_","RP$", "rp$",
           "_rp_", "\\.rp\\.", "_RF", "RIBOSEQ"),
         c("QTI"),
         c("CAGE", "cage"),
         c("80S","LSU"),
         c("40S","SSU"),
         c("ATAC"),
         c("PRPF"),
         c("PAS-Seq"),
         c("PAL-Seq"),
         c("RIP-Seq", "RIP-seq", "^rip_"),
         c("SHAPE", "Shape "),
         c("ChIP"),
         c("PAR-CLIP", "CLIP-Seq", "CLIP-seq", "iCLIP"),
         c("tRNA"),
         c("miRNA"),
         c("GRO-seq", "GROseq", "PRO-seq", "PROseq")
    )
  dt <- data.table(mainName, allNames)
  return(dt)
}

#' Get stage name variants
#'
#' Used to standardize nomeclature for experiments.\cr
#' Example: Find timepoints 2 hours, 4 hours etc.
#' Example: If using zebrafish stages as TRUE,
#' 64Cell stage is same as 2 hours post fertilization,
#' so all 2hpf will be converted to 64Cell etc.
#' @family experiment_naming
#' @param zebrafish.stages logical, FALSE. If true, convert time points to stages.
#' @references https://www.mbl.edu/zebrafish/files/2013/03/Kimmel_stagingseries1.pdf
#' @return a data.table with 2 columns, the main name, and all name variants
#' of the main name in second column as a list.
stageNames <- function(zebrafish.stages = FALSE) {
  if (zebrafish.stages) {
    mainName <- c("unfertilized", "fertilized",
                  "2to4Cell", "4Cell", "8Cell", "64Cell", "256Cell", "512Cell",
                  "1KCell", "High", "Oblong", "Sphere", "Dome", "epiboly","Shield", "Bud",
                  "Somite", "24hpf", "prim6", "prim10", "prim12",  "prim20",
                  "2dpf", "3dpf", "4dpf", "5dpf", "6dpf", "10dpf", "21dpf",
                  "24dpf")
  } else {
    mainName <- c("0h", "10min",
                  "40min", "1h", "1h15min", "2h", "2h30min", "2h45min",
                  "3h", "3h20min", "3h40min", "4h", "4h20min", "5h","6h", "10h",
                  "12h", "24h", "25h", "27h", "28h",  "33h",
                  "2d", "3d", "4d", "5d", "6d", "10d", "21d",
                  "24d")
  }

  allNames <-
    list(c("unfertilized", "Unfertilized", "_0h", "_00h", "0hpf", "^0h_"),
         c("_fertilized", "_Fertilized"),
         c("2to4Cell", "2to4cell", "2to4_cell", "2-4cell", "2-4Cell", "2-4_cell"),
         c("4cell", "4Cell", "4_cell", "_1h", "_01h", "1hpf", "^1h_"),
         c("8cell", "8Cell", "8_cell"),
         c("64cell", "64Cell", "64-cell", "64_cell", "64_Cell", "_2h", "_02h", "2hpf", "^2h_"),
         c("256cell", "256Cell", "256_cell", "256_Cell"),
         c("512cell", "512Cell", "512_cell"),
         c("1Kcell", "1KCell", "1K_cell", "_3h", "_03h", "3hpf", "^3h_"),
         c("High"),
         c("Oblong", "oblong"),
         c("Sphere", "sphere", "_4h", "_04h","4hpf", "^4h_"),
         c("Dome", "dome"),
         c("epiboly", "_5h", "_05h","5hpf", "^5_"),
         c("Shield", "shield", "_6h", "_06h", "6hpf", "^6h_"),
         c("Bud", "bud", "_10h", "10hpf", "^10h_"),
         c("Somite", "somite", "_12h", "12hpf", "^12h_"),
         c("24hpf", "_24h", "1dfp", "^24h_", " 24h "),
         c("prim6", "prim_6", "25hpf", "_25h", "^25h_"),
         c("prim10", "prim_10", "27hpf", "_27h", "^27h_"),
         c("prim12", "prim_12", "28hpf", "_28h", "^28h_"),
         c("prim20", "prim_20", "33hpf", "_33h", "^33h_"),
         c("2dpf", "_48h", "_48hpf", "^48h_", " 48h "),
         c("3dpf", " 72h "),
         "4dpf", "5dpf", "6dpf", "10dpf", "21dpf", "24dpf"
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
    "colon", "cortex", "eye", "frontal_lobe", "heart",
    "kidney", "liver", "lung", "melanocyte", "mesenchymal", "monocyte", "muscle",
    "myeloid","ovary", "prostate",
    "rectum", "retina","testis","urunary", "vagina", "skin", "tongue", "throat")
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
         c("skin", "Skin", "Dermal", "dermal"),
         c("tongue", "Tongue"),
         c("throat")
    )

  dt <- data.table(mainName, allNames)
  return(dt)
}
#' Get cell-line name variants
#'
#' Used to standardize nomeclature for experiments.\cr
#' Example: THP1 is main naming, but a variant is THP-1
#' THP-1 will then be renamed to THP1 (variables in R,
#' can not have - in them)
#' @family experiment_naming
#' @param convertToTissue logical, FALSE. If TRUE, return tissue type. NONE is
#' returned for general non-differentiated cell lines like 3T3.
#' @return a data.table with 2 columns, the main name, and all name variants
#' of the main name in second column as a list.
cellLineNames <- function(convertToTissue = FALSE) {
  if (convertToTissue) {
    mainName <- c("ovary", "lung", "lung","breast", "kidney", "ovary",
                  "liver", "NONE", "colon", "breast", "blood", "liver",
                  "NONE", "bone","breast", "breast", "blood","throat",
                  "pancreas", "blood", "skin", "breast", "bone", "prostate",
                  "NONE")
  } else {
    mainName <- c("A2780", "A549", "Calu3","CN34", "HEK293", "HeLa",
                  "HepG2", "Hesc", "HCT116", "Hmel", "HSB2", "Huh7",
                  "iPSC", "K562","MCF10A", "MDA", "MV4-11","NPC",
                  "PANC1","THP1", "TSC2", "T47D","U2OS", "PC3",
                  "3T3")
  }


  allNames <-
    list(c("A2780"),
         c("A549"),
         c("Calu 3", "Calu3"),
         c("CN34"),
         c("HEK293", "Hek293", "293T", "HEK_", "T-REx"),
         c("HeLa", "HELA", "Hela", "hela"),
         c("HepG2"),
         c("hesc", "Hesc"),
         c("HCT116"),
         c("Hmel"),
         c("HSB2"),
         c("Huh7"),
         c("iPSC", "ipsc"),
         c("K562"),
         c("MCF10-A", "MCF10A"),
         c("MDA-"),
         c("MV4-11"),
         c("npc", "NPC"),
         c("PANC1"),
         c("THP-1", "THP1"),
         c("TSC2"),
         c("T47D"),
         c("U2OS"),
         c("PC3", "PC-3"),
         c("3T3")
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
    list(c("rep1", " rep 1 ", "Rep1", "rep-1","replicate1", "run1", "run_1_", "_r1_", "WT1", " 1$", "_01$"),
         c("rep2", " rep 2 ", "Rep2", "rep-2","replicate2", "run2", "run_2_", "_r2_", "WT2", " 2$", "_02$"),
         c("rep3", " rep 3 ", "Rep3", "rep-3","replicate3", "run3", "run_3_", "_r3_", "WT3", " 3$", "_03$"),
         c("rep4", " rep 4 ", "Rep4", "rep-4","replicate4", "run4", "run_4_", "_r4_", "WT4", " 4$", "_04$"),
         c("rep5", " rep 5 ", "Rep5", "rep-5","replicate5", "run5", "run_5_", "_r5_", "WT5", " 5$", "_05$"),
         c("rep6", " rep 6 ", "Rep6", "rep-6","replicate6", "run6", "run_6_", "_r6_", "WT6", " 6$", "_06$")
    )
  dt <- data.table(mainName, allNames)
  return(dt)
}

#' Get batch name variants
#'
#' Used to standardize nomeclature for experiments.\cr
#' Example: Biological samples (batches) batch will become b1
#' @family experiment_naming
#' @return a data.table with 2 columns, the main name, and all name variants
#' of the main name in second column as a list.
batchNames <- function() {
  mainName <- c("b1", "b2", "b3", "b4", "b5", "b6")
  allNames <-
    list(c("batch1", "batch 1 ", "batch 1)", "set1"),
         c("batch2", "batch 2 ", "batch 2)", "set2"),
         c("batch3", "batch 3 ", "batch 3)", "set3"),
         c("batch4", "batch 4 ", "batch 4)", "set4"),
         c("batch5", "batch 5 ", "batch 5)", "set5"),
         c("batch6", "batch 6 ", "batch 6)", "set6")
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
  mainName <- c("WT", "MZ", "4Ei",
                "Mutant", "cas9", "NMDA", "DHPG", "KD",
                "KO", "DKO", "high", "low")
  allNames <-
    list(c("WT", "wt", "wild_type", "Wild_type", "wild-type",
           "untreated", "control", "Control", "CTRL", "Ctrl","Basal", "_mock_",
           " mock "),
         c("MZ", "dicer"),
         c("4Ei", "4ei"),
         c("mutant", "Mutant"),
         c("Cas9", "cas9"),
         c("NMDA"),
         c("DHPG"),
         c("knockdown", "Knockdown", "_kd_", " KD ", "-KD ", "_KD_"),
         c("knockout", "Knockout", " KO ", "-KO ", "_KO_", "-KO-"),
         c("double knockout", " DKO ", "_DKO_"),
         c("high", "\\+"),
         c("low")
    )
  dt <- data.table(mainName, allNames)
  return(dt)
}

#' Get translocation inhibitor name variants
#'
#' Used to standardize nomeclature for experiments.\cr
#' Example: cycloheximide, lactimidomycin, harringtonine
#' @family experiment_naming
#' @return a data.table with 2 columns, the main name, and all name variants
#' of the main name in second column as a list.
inhibitorNames <- function() {
  mainName <- c("chloram","chx", "harr", "frozen", "lactim")
  allNames <-
    list(c("chloramphenicol", "Chloramphenicol"),
         c("cycloheximide", "Cycloheximide", "chx", "CHX"),
         c("harringtonin", "Harringtonin", "_har_", "_harr_", "_harr$", "_Harr$"),
         c("frozen", "freeze", "freezing"),
         c("lactimidomycin", "Lactimidomycin", "ltm", "_LTM$")
    )
  dt <- data.table(mainName, allNames)
  return(dt)
}

#' Get cell fraction name variants
#'
#' Used to standardize nomeclature for experiments.\cr
#' Example: cytosolic, mitochondrial, specific gene knock down
#' @family experiment_naming
#' @return a data.table with 2 columns, the main name, and all name variants
#' of the main name in second column as a list.
fractionNames <- function() {
  mainName <- c("cyto","mito", "nuc","dmso", "thaps", "auxin", "silvestrol")
  allNames <-
    list(c("cyto"),
         c("mito"),
         c("_nuc_", "nuclear"),
         c("DMSO", "dmso"),
         c("Thapsigargin", "thapsigargin"),
         c("auxin"),
         c("silvestrol")
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
