#!/bin/bash

#HT 12/02/19
# script to process and align genomic fastq datasets
#1 STAR INDICES (If not existing, must be done in seperate script for now)
#2 make directories
#3 Trim adaptor
#4 Remove PhiX
#5 Remove rRNA (from silva)
#6 Remove organism specific ncRNA (from ensembl, zebrafish does not contain tRNA)
#7 Remove organism specific tRNA (from tRNAscan-SE)
#8 Map to organism specific reference

usage(){
cat << EOF
usage: $0 options

script to process and align genomic single end or paired end fastq datasets

OPTIONS:
	Important options:
	-f	path to input fasta file. Also define -F if paired end! Must be file type of: fasta, fa, fastq, fq or gz
	-F 	path to input fastq file 2 (paired end reads 2) Must be file type of: fasta, fa, fastq, fq or gz
	-o	path to output dir
	-l	minimum length of reads (default: 20)
	-T	max mismatches of reads (default: 3)
	-g	genome dir for all indices (Standard is zebrafish: danrerio10, change to human index if needed etc)
  -s	steps of depletion and alignment wanted:
		(a string: which steps to do? (default: "tr-ge", write "all" to get all: "tr-ph-rR-nc-tR-ge",
		   or tr-co-ge, depending on if you merged contaminants or not.)
			 tr: trim, co: contaminants, ph: phix, rR: rrna, nc: ncrna, tR: trna, ge: genome)
		Write your wanted steps, seperated by "-". Order does not matter.
		To just do trim and alignment to genome write -s "tr-ge"
	-a	adapter sequence for trim (found automaticly if not given), also you can write -a "disable",
		to disable it or "standard" to get "AAAAAAAAAA", the illumina standard sequence.
	-t	trim front (default 3) How many bases to pre trim reads 5' end,
	        as it frequently represents an untemplated addition during reverse transcription.
	-A	Alignment type: (default Local, EndToEnd (Local is Local, EndToEnd is force Global))
  -B Allow introns (default yes (1), else no (0))
	Path arguments:
	-S      path to STAR (default: ~/bin/STAR-2.7.0c/source/STAR)
	-P      path to fastp (trimmer) (default: ~/bin/fastp)

	Less important options:
	-r	resume?: a character (defualt n) (n for new start fresh with file f from point s,
			             (if you want a continue from crash specify the step you want to start
				      from, as defined in -s, start on genome, do -r "ge")
	-m	max cpus allowed (defualt 90)
	-M  Max multimapping (default 10) Set to 1 to get only unique reads. Only applies for genome
      step, not the depletion step.
	-k	a character, Keep loaded genomes STAR index:
	    yes (y), no: Remove loaded (n), no shared genome (noShared),
	    default (n)
	-K  Keep contaminant aligned bam files. Default: "no", alternative: "yes".
	-q	a character, Do quality filtering:
	    yes: "default" no: "disable". Uses default fastp QF.
	-u  Keep unaligned reads from genome alignment step.
	      Default "None" (no), alternative: "Fastx""
	-h	this help message

fastp location must be: ~/bin/fastp
(if you don't have it install to bin folder from: https://github.com/OpenGene/fastp

(if you don't have it install to bin folder from: https://github.com/alexdobin/STAR

NOTES:
if STAR is stuck on load, run this line:
STAR --genomeDir /references/human/STAR_INDEX/genomeDir/ --genomeLoad Remove

example usage: RNA_Align_pipeline.sh -f <in.fastq.gz> -o <out_dir>

EOF
}

# Default arguments:
min_length=20
mismatches=3
gen_dir=""
allSteps="tr-ge"
steps=$allSteps
resume="n"
alignment="Local"
allow_introns=0
adapter="auto"
quality_filtering="disable"
maxCPU=90
multimap=10
trim_front=3
keep="n"
keepContam="no"
keep_unmapped_genome="None"
STAR="~/bin/STAR-2.7.0c/source/STAR"
fastp="~/bin/fastp"
while getopts ":f:F:o:l:T:g:s:a:t:A:B:r:m:M:K:k:p:S:P:q:u:h" opt; do
    case $opt in
    f)
        in_file=$OPTARG
        echo "-f input file: $OPTARG"
	      ;;
    F)
      	in_file_two=$OPTARG
      	echo "-F input file 2: $OPTARG"
	      ;;
    o)
        out_dir=$OPTARG
        echo "-o output folder: $OPTARG"
        ;;
    l)
        min_length=$OPTARG
        echo "-l minimum length of reads: $OPTARG"
        ;;
    T)
        mismatches=$OPTARG
        echo "-T max mismatches of reads: $OPTARG"
        ;;
    g)
        gen_dir=$OPTARG
        echo "-g genome dir for all indices: $OPTARG"
        ;;
    s)
        steps=$OPTARG
        echo "-s steps to do: $OPTARG"
        ;;
    a)
        adapter=$OPTARG
        echo "-a adapter sequence: $OPTARG"
        ;;
    q)
	      quality_filtering=$OPTARG
	      echo "-q quality filtering: $quality_filtering"
        ;;
    t)
        trim_front=$OPTARG
        echo "-t trim front (nt): $OPTARG"
        ;;
    A)
        alignment=$OPTARG
        echo "-A alignment type: $OPTARG"
        ;;
    B)
        allow_introns=$OPTARG
        echo "-B allow_introns: $OPTARG"
        ;;
    r)
	      resume=$OPTARG
	      echo "-r resume (r or new n): $OPTARG"
        ;;
    m)
      	maxCPU=$OPTARG
      	echo "-m maxCPU: $OPTARG"
        ;;
    M)
      	multimap=$OPTARG
      	echo "-M max multimap: $OPTARG"
        ;;
    S)
      	STAR=$OPTARG
      	echo "-S STAR location: $OPTARG"
        ;;
    P)
      	fastp=$OPTARG
      	echo "-P fastp location: $OPTARG"
        ;;
    k)
      	keep=$OPTARG
      	echo "-k Keep Star Index loaded: $OPTARG"
        ;;
    K)
      	keepContam=$OPTARG
      	echo "-K Keep contamination reads: $OPTARG"
        ;;
    u)
      	keep_unmapped_genome=$OPTARG
      	echo "-u Keep unmapped genome reads: $OPTARG"
        ;;
    h)
        usage
        exit
        ;;
    ?)
        echo "Invalid option: -$OPTARG"
        usage
        exit 1
        ;;
    esac
done
echo ""
#exit 1

if [ "$steps" == "all" ]; then
	steps="tr-ph-rR-nc-tR-ge"
fi

IFS="-" read -a stepsArray <<< "$steps"
export stepsArray


if [ -z "$out_dir" ]; then
	echo "Error, out directory (-o) must be speficied!"
	exit 1
fi

if [ ! -z "$in_file" ]; then
	if [ ! -f "$in_file" ]; then
	    echo "input file f does not name existing file!"
	    exit 1
	fi

	if [[ ! "$in_file" =~ .*\.(fasta|fa|fastq|gz|fq) ]]; then
	    echo "Invalid input file type: $in_file"
	    echo "Must be either fasta, fa, fastq, fq or gz"
	    exit 1
	fi
else
	echo "Error input file :f was not assigned!"
	exit 1
fi

if [ ! -z "$in_file_two" ]; then

	if [ ! -f "$in_file_two" ]; then
	    echo "Error input file 2 :F does not name existing file!"
	    exit 1
	fi

	if [[ ! $in_file_two =~ .*\.(fasta|fa|fastq|gz|fq) ]]; then
	    echo "Invalid input file type: $in_file_two"
	    echo "Must be either fasta, fa, fastq, fq or gz"
	    exit 1
	fi
fi


# 1. mkdir

if [ ! -d "$out_dir" ]; then
    mkdir -p $out_dir
fi

contamSAMmode="None"
contamSAMtype="None"
if [ $keepContam == "yes" ]; then
  contamSAMmode="Full"
  contamSAMtype="BAM Unsorted"
fi

# 3-7 filtering
contaminants=$gen_dir/contaminants_genomeDir
phix=$gen_dir/PhiX_genomeDir
rRNA=$gen_dir/rRNA_genomeDir
ncRNA=$gen_dir/ncRNA_genomeDir
tRNA=$gen_dir/tRNA_genomeDir
usedGenome=$gen_dir/genomeDir
if [ ! -d "$usedGenome" ]; then
      if [ $steps == "tr" ]; then
        echo "Running trim only mode"
      else
		    echo "Error: the given STAR index dir does not exist!"
	      exit 1
	    fi
fi

ibn=$(basename ${in_file}) # <- name to use
ibn=${ibn%.gz}
ibn=${ibn%.fastq}
ibn=${ibn%.fq}
ibn=${ibn%.fasta}
ibn=${ibn%.fa}
echo "basename of file is: ${ibn}"

# Index of current step (tr-ge, means tr = 0)
function indexOf()
{
	string=($1) && shift
	myArray=($@)

	for i in "${!myArray[@]}"; do

		if [ ${myArray[i]} == $string ]; then
        		echo $i
		fi
    	done
}

# Get all relative paths to outputs of pipeline
# 1 ${stepsArray[ind]} 2 ${out_dir} 3 ${ibn} 4 ${in_file_two}
# Return: string Full path to designated file (2 for paired end)
function pathList()
{
	if [ -z ${4} ]; then # if single end
		case $1 in
		  "co")
		    echo "${2}/contaminants_depletion/contaminants_${3}_Unmapped.out.mate1"
		    ;;
			"tR")
				echo "${2}/tRNA_depletion/tRNA_${3}_Unmapped.out.mate1"
				;;
			"nc")
				echo "${2}/ncRNA_depletion/ncRNA_${3}_Unmapped.out.mate1"
				;;
			"rR")
				echo "${2}/rRNA_depletion/rRNA_${3}_Unmapped.out.mate1"
				;;
			"ph")
				echo "${2}/phix_depletion/PhiX_${3}_Unmapped.out.mate1"
				;;
			"tr")
				echo "${2}/trim/trimmed_${3}.fastq"
				;;
		esac
	else # if paired end
		case $1 in
		  "co")
				echo "${2}/contaminants_depletion/contaminants_${3}_Unmapped.out.mate1 ${2}/contaminants_depletion/contaminants_${3}_Unmapped.out.mate2"
				;;
			"tR")
				echo "${2}/tRNA_depletion/tRNA_${3}_Unmapped.out.mate1 ${2}/tRNA_depletion/tRNA_${3}_Unmapped.out.mate2"
				;;
			"nc")
				echo "${2}/ncRNA_depletion/ncRNA_${3}_Unmapped.out.mate1 ${2}/ncRNA_depletion/ncRNA_${3}_Unmapped.out.mate2"
				;;
			"rR")
				echo "${2}/rRNA_depletion/rRNA_${3}_Unmapped.out.mate1 ${2}/rRNA_depletion/rRNA_${3}_Unmapped.out.mate2"
				;;
			"ph")
				echo "${2}/phix_depletion/PhiX_${3}_Unmapped.out.mate1  ${2}/phix_depletion/PhiX_${3}_Unmapped.out.mate2"
				;;
			"tr")
				echo "${2}/trim/trimmed_${3}.fastq ${2}/trim/trimmed2_${3}.fastq"
				;;
		esac
	fi
}

# Get fasta/bam file to use in this step
# Parameters 1. $resume 2. $in_file 3.current:'ge'
# 4. ${out_dir} 5. ${ibn} 6. ${in_file_two}
# Return: string: 1 path for single end, 2 for paired.
function inputFile()
{
	var=$(indexOf "$3" ${stepsArray[@]})
	if [[ "$var" == "0" ]]; then
		echo "${2} ${6}"
	else
		ind=$(expr $var - 1)
		echo $(pathList ${stepsArray[ind]} $4 $5 ${6})
	fi
}

# Decide to do this step or not, given input
# Parameters $1 resume (y, n), $2 current:'ge', $3 $steps tr-gr
# TODO: fix for $1 == "c", the 1st order else statement
function doThisStep()
{
  if [ $1 == "n" ]; then
  	if grep -q $2 <<< $3; then
  		echo "yes"
  	else
  		echo "no"
  	fi
  else
  	if [ $2 == $1 ]; then
  		echo "yes"
  	else
  		echo "no"
  	fi
  fi
}

# Should STAR use zcat or no decompression ?
# Return: string zcat or "-"
function comp()
{
	if [[ $1 =~ .*\.(gz) ]]; then
		echo "zcat"
	else
		echo "-"
	fi
}

# 1: max cores, 2: currently used cores
function nCores()
{
	if (( $1 > $2 )); then
		echo $2
	else
		echo $1
	fi
}

# Keep loaded STAR index of genomes (y) or not (n), default (n)
function keepOrNot()
{
	if [[ "$1" == "y" ]]; then
		echo LoadAndKeep
	elif [[ "$1" == "n" ]]; then
		echo LoadAndRemove
	elif [[ "$1" == "noShared" ]]; then
	  echo NoSharedMemory
	else
		exit 1
	fi
}

# Add paired read 2, if exists
#  ${out_dir} ${ibn} ${in_file_two}
function trimPaired()
{
	if [ ! -z "$3" ]; then
		echo "${1}/trim/trimmed2_${2}.fastq"
	else
		echo ""
	fi
}

# Adapter definitions for fastp
if [ -z  "$adapter" ]; then
	adapter="" # "" is auto detection
elif [ $adapter == "auto" ]; then
	adapter="" # "auto" is auto detection
elif [ $adapter == "autoPE" ]; then
	adapter="--detect_adapter_for_pe" # auto detection for PE
elif [ $adapter == "disable" ]; then
	adapter="--disable_adapter_trimming"
else
	# Check if it is one of the templates
	if [ $adapter == "standard" ]; then
		adapter="AAAAAAAAAA"
	elif [ $adapter == "illumina" ]; then
	  echo "Using Illumina preset adapter"
		adapter="AGATCGGAAGAGC"
	elif [ $adapter == "small_RNA" ]; then
	  echo "Using Small RNA preset adapter"
		adapter="TGGAATTCTCGG"
	elif [ $adapter == "nextera" ]; then
	  echo "Using Nextera preset adapter"
		adapter="CTGTCTCTTATA"
	elif [ $adapter == "ingolia12" ]; then
	  echo "Using Ingolia(2012 Ribo-seq) preset adapter"
		adapter="CTGTAGGCACCATCAAT"
	fi
	adapter="--adapter_sequence=${adapter}"
fi

# Quality filtering for fastp
if [ $quality_filtering == "default" ]; then
	quality_filtering="" # Default QF
elif [ $quality_filtering == "disable" ]; then
	quality_filtering="--disable_quality_filtering"
fi

#------------------------------------------------------------------------------------------
    #3 FASTP (Trim adaptors, Cut, Quality, fastq report)
    #------------------------------------------------------------------------------------------
    #--in1 input file
    #--out1 output name
    #--trim_front1 trim 3 bases
    #--length_required minimum length
    #--disable_quality_filtering no fastq filtering
    #--adapter_sequence adapter sequence, normally set manually (normally needed)
if [ $(doThisStep $resume 'tr' $steps) == "yes" ]; then
	echo trimming
	if [ ! -d ${out_dir}/trim ]; then
        mkdir ${out_dir}/trim
        if [ ! -d ${out_dir}/trim ]; then
          echo "Error: could not create trim dir, do you have access to disc?"
          exit 1
        fi
  fi

	eval $fastp \
		--in1=${in_file} \
		--in2="${in_file_two}" \
		--out1=${out_dir}/trim/trimmed_${ibn}.fastq \
		--out2=$(trimPaired ${out_dir} ${ibn} ${in_file_two}) \
		--json=${out_dir}/trim/report_${ibn}.json \
		--html=${out_dir}/trim/report_${ibn}.html \
		--trim_front1=${trim_front} \
		--trim_front2=${trim_front} \
		--length_required=$min_length \
		$quality_filtering \
		$adapter \
		--thread $(nCores 16 $maxCPU)
fi

#------------------------------------------------------------------------------------------
    #4 (alternative): Remove merged contaminants
    #------------------------------------------------------------------------------------------
# get output of everything that did not hit, as fastq
if [ $(doThisStep $resume 'co' $steps) == "yes" ]; then
	echo "Contaminant depletion:"
	if [ ! -d ${out_dir}/contaminants_depletion ]; then
        mkdir ${out_dir}/contaminants_depletion
  fi

	eval $STAR \
	--readFilesIn $(inputFile $resume $in_file "co" ${out_dir} ${ibn} ${in_file_two}) \
	--genomeDir ${contaminants} \
	--genomeLoad $(keepOrNot $keep)  \
	--outFileNamePrefix ${out_dir}/contaminants_depletion/contaminants_${ibn}_ \
	--outSAMtype $contamSAMtype \
	--outSAMmode $contamSAMmode \
	--outReadsUnmapped Fastx \
	--outFilterMatchNmin $min_length \
	--runThreadN $(nCores 90 $maxCPU) \
	--readFilesCommand $(comp $(inputFile $resume $in_file 'co' ${out_dir} ${ibn})) \
	--limitIObufferSize 50000000
fi

 #------------------------------------------------------------------------------------------
    #4 Remove PhiX
    #------------------------------------------------------------------------------------------
    #--readFilesIn input file
    #--genomeDir STAR index genome dir
    #--outFileNamePrefix output prefix name (will add input name also)
    #--outSAMtype type of SAM output (bam, sam or None)
    #--outReadsUnmapped output type for unmapped reads
    #--outFilterMatchNmin minimum length of reads accepted
    #--runThreadN number of threads to use
    #--limitIObufferSize hard drive buffer size (smaller is better when IO is bottle neck)
# get output of everything that did not hit, as fastq
if [ $(doThisStep $resume 'ph' $steps) == "yes" ]; then
	echo "PhiX depletion:"
	if [ ! -d ${out_dir}/phix_depletion ]; then
        mkdir ${out_dir}/phix_depletion
  fi

	eval $STAR \
	--readFilesIn $(inputFile $resume $in_file "ph" ${out_dir} ${ibn} ${in_file_two}) \
	--genomeDir ${phix} \
	--genomeLoad $(keepOrNot $keep) \
	--outFileNamePrefix ${out_dir}/phix_depletion/PhiX_${ibn}_ \
	--outSAMtype $contamSAMtype \
	--outSAMmode $contamSAMmode \
	--outReadsUnmapped Fastx \
	--outFilterMatchNmin $min_length \
	--runThreadN $(nCores 70 $maxCPU) \
	--readFilesCommand $(comp $(inputFile $resume $in_file 'ph' ${out_dir} ${ibn})) \
	--limitIObufferSize 50000000 \
	--alignIntronMax 1
fi


#------------------------------------------------------------------------------------------
    #5 Remove rRNA
    #------------------------------------------------------------------------------------------
# get output of everything that did not hit, as fastq
if [ $(doThisStep $resume 'rR' $steps) == "yes" ]; then
	echo "rRNA depletion:"
	if [ ! -d ${out_dir}/rRNA_depletion ]; then
        mkdir ${out_dir}/rRNA_depletion
  fi

	eval $STAR \
	--readFilesIn $(inputFile $resume $in_file "rR" ${out_dir} ${ibn} ${in_file_two}) \
	--genomeDir ${rRNA} \
	--genomeLoad $(keepOrNot $keep)  \
	--outFileNamePrefix ${out_dir}/rRNA_depletion/rRNA_${ibn}_ \
	--outSAMtype $contamSAMtype \
	--outSAMmode $contamSAMmode \
	--outReadsUnmapped Fastx \
	--outFilterMatchNmin $min_length \
	--runThreadN $(nCores 90 $maxCPU) \
	--readFilesCommand $(comp $(inputFile $resume $in_file 'rR' ${out_dir} ${ibn})) \
	--limitIObufferSize 50000000 \
	--alignIntronMax 1
fi

 #------------------------------------------------------------------------------------------
    #6 Remove organism specific ncRNA
    #------------------------------------------------------------------------------------------
# get output of everything that did not hit, as fastq
if [ $(doThisStep $resume 'nc' $steps) == "yes" ]; then
	echo "ncRNA depletion:"
	if [ ! -d ${out_dir}/ncRNA_depletion ]; then
    	mkdir ${out_dir}/ncRNA_depletion
  fi

	eval $STAR \
	--readFilesIn $(inputFile $resume $in_file "nc" ${out_dir} ${ibn} ${in_file_two})\
	--genomeDir ${ncRNA} \
	--genomeLoad $(keepOrNot $keep)  \
	--outFileNamePrefix ${out_dir}/ncRNA_depletion/ncRNA_${ibn}_ \
	--outSAMtype $contamSAMtype \
	--outSAMmode $contamSAMmode \
	--outReadsUnmapped Fastx \
	--outFilterMatchNmin $min_length \
	--runThreadN $(nCores 80 $maxCPU) \
	--readFilesCommand $(comp $(inputFile $resume $in_file 'nc' ${out_dir} ${ibn})) \
	--limitIObufferSize 50000000
fi

#------------------------------------------------------------------------------------------
    #7 Remove organism specific tRNA
    #------------------------------------------------------------------------------------------
# get output of everything that did not hit, as fastq
# ENCODE tRNA mapping defaults
if [ $(doThisStep $resume 'tR' $steps) == "yes" ]; then
	echo "tRNA depletion"
	if [ ! -d ${out_dir}/tRNA_depletion ]; then
        mkdir ${out_dir}/tRNA_depletion
  fi

	eval $STAR \
	--readFilesIn $(inputFile $resume $in_file "tR" ${out_dir} ${ibn} ${in_file_two}) \
	--genomeDir ${tRNA} \
	--genomeLoad $(keepOrNot $keep)  \
	--outFileNamePrefix ${out_dir}/tRNA_depletion/tRNA_${ibn}_ \
	--outSAMtype $contamSAMtype \
	--outSAMmode $contamSAMmode \
	--outReadsUnmapped Fastx \
	--outFilterMatchNmin $min_length \
	--runThreadN $(nCores 70 $maxCPU) \
  --readFilesCommand $(comp $(inputFile $resume $in_file 'tr' ${out_dir} ${ibn})) \
	--limitIObufferSize 50000000 \
	--alignIntronMax 1 \
	--seedPerWindowNmax 20 \
	--outFilterMultimapNmax 20
fi
#------------------------------------------------------------------------------------------
    # 8. aligner (STAR)
    #------------------------------------------------------------------------------------------
    #--limitBAMsortRAM RAM used by sorting function (higher is better)
    #--alignEndsType how to align (local alignment or whole read(harder if adapter or weak 3' end))
    # <(gunzip -c ${in_file})
if [ $(doThisStep $resume 'ge' $steps) == "yes" ]; then
	echo "Final mapping to genome:"
	if [ ! -d ${out_dir}/aligned ]; then
        mkdir ${out_dir}/aligned
  fi
  ((allow_introns ^= 1)) # XOR to flip, since STAR 0 is allow introns

	eval $STAR \
	--readFilesIn $(inputFile $resume $in_file 'ge' ${out_dir} ${ibn} ${in_file_two}) \
	--genomeDir ${usedGenome} \
	--outFileNamePrefix ${out_dir}/aligned/${ibn}_ \
	--outSAMtype BAM SortedByCoordinate \
	--outReadsUnmapped $keep_unmapped_genome \
	--runThreadN $(nCores 80 $maxCPU) \
	--genomeLoad $(keepOrNot $keep)  \
	--limitIObufferSize 50000000 \
	--outFilterMatchNmin $min_length \
	--limitBAMsortRAM 30000000000 \
	--readFilesCommand $(comp $(inputFile $resume $in_file 'ge' ${out_dir} ${ibn})) \
	--alignEndsType $alignment \
	--alignIntronMax $allow_introns \
	--outFilterMultimapNmax $multimap \
	--outFilterMismatchNmax $mismatches
fi

#TODO
# Remove empty folders, and possibly logs in seperate folder ?
