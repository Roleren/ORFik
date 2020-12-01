#!/bin/bash

#HT 12/02/19
# script to process and align genomic non paired fastq datasets
#1 STAR INDICES (If not existing, must be done in seperate script for now)
#2 make directories
#3 Trim adaptor
#4 (alternative:) Remove merged contaminants
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
	-f	path to input folder (only fasta/q files in folder allowed for now)
		Must be file types of: fasta, fa, fastq, fq or gz (single or paired end reads)
	-o	path to output dir
	-p	paired end? (yes, defualt: no)
	-l	minimum length of reads (default: 20)
	-T	max mismatches of reads (default: 3)
	-g	genome dir for all STAR indices
	-s	steps of depletion and alignment wanted:
		(a string: which steps to do? (default: "tr-ge", write "all" to get all: "tr-ph-rR-nc-tR-ge",
		   or tr-co-ge, depending on if you merged contaminants or not.)
			 tr: trim, co: contaminants, ph: phix, rR: rrna, nc: ncrna, tR: trna, ge: genome)
		Write your wanted steps, seperated by "-". Order does not matter.
		To just do trim and alignment to genome write -s "tr-ge"
	-a	adapter sequence for trim (found automaticly if not given), also you can write -a "disable",
		to disable it
	-t	trim front (default 3) How many bases to pre trim reads 5' end,
	        as it frequently represents an untemplated addition during reverse transcription.
	-A	Alignment type: (default Local, EndToEnd (Local is Local, EndToEnd is force Global))
  -M  Max multimapping (default 10) Set to 1 to get only unique reads. Only applies for genome
      step, not the depletion step.
	Path arguments:
	-S   path to STAR (default: ~/bin/STAR-2.7.0c/source/STAR)
	-P   path to fastp (trimmer) (default: ~/bin/fastp)
	-C	 path to cleaning script, internal

	Less important options:
	-r	resume?: a character (defualt n) (n for new start fresh with file f from point s,
			             (if you want a continue from crash specify the step you want to start
				      from, as defined in -s, start on genome, do -r "ge")
	-m	max cpus allowed (defualt 90)
	-i	include subfolders (defualt n, for no), if you want subfolder do y, for yes.
	-h	this help message

fastp location must be: ~/bin/fastp
STAR location must be: ~/bin/STAR-2.7.0c/source/STAR

NOTE: if STAR is stuck on load, run this line:
STAR --genomeDir /export/valenfs/data/references/Zv10_zebrafish_allsteps --genomeLoad Remove

example usage: RNA_Align_pipeline_folder.sh -f <in.fastq.gz> -o <out_dir>

EOF
}

# Default arguments:
echo "##############################################"
echo $'\nArguments for folder run are the following:'
min_length=20
mismatches=3
gen_dir=""
allSteps="tr-ge"
steps=$allSteps
resume="n"
alignment="Local"
adapter="auto"
maxCPU=90
multimap=10
subfolders="n"
trim_front=3
paired="no"
STAR="~/bin/STAR-2.7.0c/source/STAR"
fastp="~/bin/fastp"
align_single=""
cleaning=""
while getopts ":f:o:p:l:T:g:s:a:t:A:r:m:M:S:i:P:I:C:h:" opt; do
    case $opt in
    f)
        in_dir=$OPTARG
        echo "-f input folder $OPTARG"
	      ;;
    o)
        out_dir=$OPTARG
        echo "-o output folder $OPTARG"
        ;;
    p)
	      paired=$OPTARG
        echo "-p paired end $OPTARG"
        ;;
    l)
        min_length=$OPTARG
        echo "-l minimum length of reads $OPTARG"
        ;;
    T)
        mismatches=$OPTARG
        echo "-T max mismatches of reads $OPTARG"
        ;;
    g)
        gen_dir=$OPTARG
        echo "-g genome dir for all STAR indices $OPTARG"
        ;;
    s)
        steps=$OPTARG
        echo "-s steps to do: $OPTARG"
        ;;
    a)
        adapter=$OPTARG
        echo "-a adapter sequence $OPTARG"
        ;;
    t)
        trim_front=$OPTARG
        echo "-t trim front number $OPTARG"
        ;;
    A)
        alignment=$OPTARG
        echo "-A alignment type $OPTARG"
        ;;
    r)
	resume=$OPTARG
	echo "-r resume (new: n, or step as ge): $OPTARG"
        ;;
    m)
	maxCPU=$OPTARG
	echo "-m maxCPU: $OPTARG"
        ;;
    M)
	multimap=$OPTARG
	echo "-m max multimap: $OPTARG"
        ;;
    i)
	subfolders=$OPTARG
	echo "-i subfolders: $OPTARG"
        ;;
    S)
	STAR=$OPTARG
	echo "-S STAR location: $OPTARG"
        ;;
    P)
	fastp=$OPTARG
	echo "-P fastp location: $OPTARG"
        ;;
    C)
	cleaning=$OPTARG
	echo "-C cleaning location: $OPTARG"
        ;;
    I)
	align_single=$OPTARG
	echo "-I align_single location: $OPTARG"
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

echo $'\n'

if [ $steps == "all" ]; then
	steps="tr-ph-rR-nc-tR-ge"
fi
if [ -z $out_dir ]; then
	echo "Error, out directory (-o) must be speficied!"
	exit 1
fi

# Normal paired end read run
# 1 check if there exists a grouping, and the are equal number in each grouping
# 2 make a list for each end
# send each end into the paired end function
#TODO: Fix checks, incase files are not in order
function findPairs()
{
	f=($1) && shift
	myArray=($@)
	keep="y"
	i=0
	echo "${#myArray[@]}"
	echo $(( ${#myArray[@]} % 2 ))

	if [ $(( ${#myArray[@]} % 2 )) == 1 ]; then
		echo "folder must have even number of fasta/q files, for paired end run!"
		exit 1
	fi

	suffix="_001.fastq.gz"
        # remove suffix to get matched pairs
	#for i in "${!myArray[@]}"; do

		#bn=${myArray[i]} | sed  -e "s/$suffix$//"
		#echo ${myArray[i]} | sed  -e "s/$suffix$//"
        	#echo $(basename ${bn})

    	#done
  echo "#############################################"
	for ((x=0; x<${#myArray[@]}; x = x + 2));
	do
    		echo "running paired end for files: $f/${myArray[x]} and $f/${myArray[x+1]}"
    		echo "files  $((x + 1)) and $((x + 2)) / $numOfFiles"
		a="$f/${myArray[x]}"
		b="$f/${myArray[x+1]}"
		i=$((x + 2))
		if [[ $i == $numOfFiles ]];then
			keep="n"
		fi

		eval $align_single -o "$out_dir" -f "$a" -F "$b"  -a "$adapter" -s "$steps" -r "$resume" -l "$min_length" -T $mismatches -g "$gen_dir" -m "$maxCPU" -A "$alignment" -t "$trim_front" -k $keep -P "$fastp" -S "$STAR"
    echo "-------------------------------------------"
	done
}
#TODO: find a way to remove this
function findPairsSub()
{
	myArray=($@)
	echo "${#myArray[@]}"
	echo $(( ${#myArray[@]} % 2 ))

	if [ $(( ${#myArray[@]} % 2 )) == 1 ]; then
		echo "folder must have even number of fasta/q files, for paired end run!"
		exit 1
	fi
  echo "#############################################"
	for ((x=0; x<${#myArray[@]}; x = x + 2));
	do
    echo "running paired end with subfolders for files:\n ${myArray[x]} and\n ${myArray[x+1]}"
		a="${myArray[x]}"
		b="${myArray[x+1]}"

		eval $align_single -o "$out_dir" -f "$a" -F "$b"  -a "$adapter" -s "$steps" -r "$resume" -l "$min_length" -T $mismatches -g "$gen_dir" -m "$maxCPU" -M "$multimap" -A "$alignment" -t "$trim_front" -k "y" -P "$fastp" -S "$STAR"
		echo "-------------------------------------------"

	done
}
# Run per file / pair
listOfFiles=$(ls ${in_dir} | grep '\.fasta\|\.fa\|\.fastq\|\.fq')
numOfFiles=$(ls ${in_dir} | grep '\.fasta\|\.fa\|\.fastq\|\.fq' | wc -l)
echo "The total number of files are:"
echo $numOfFiles
i=0
keep="y"
if [ $paired == "yes" ]; then
	if [ $subfolders == "n" ]; then
		findPairs $in_dir $listOfFiles
	else
		findPairsSub $(find ${in_dir} | grep '\.fasta\|\.fa\|\.fastq\|\.fq\|\.gz' | sort)
	fi

else
	for x in $listOfFiles
	do
		echo "running single end for file: $x"
		i=$((i + 1))
		echo "file  $i / $numOfFiles"

		if [[ $i == $numOfFiles ]];then
			keep="n"
		fi
    x=$in_dir/$x

		eval $align_single -o "$out_dir" -f "$x"  -a "$adapter" -s "$steps" -r "$resume" -l "$min_length" -T $mismatches -g "$gen_dir" -m "$maxCPU" -M "$multimap" -A "$alignment" -t "$trim_front" -k $keep -P "$fastp" -S "$STAR"
		echo "----------------------------------------------"
	done
fi

### Cleanup
# Folder cleanup
eval $cleaning $out_dir
# Log command run
echo "./RNA_Align_pipeline_folder.sh $@" > $out_dir/runCommand.log
