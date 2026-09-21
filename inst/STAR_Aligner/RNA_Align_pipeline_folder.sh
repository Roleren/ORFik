#!/bin/bash
# Stop immediately on failed processing commands, including pipeline failures.
set -eo pipefail

#HT 12/02/19
# script to process and align genomic fastq datasets
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
		Write your wanted steps, seperated by "-". Use processing order, with trimming first and genome alignment last.
		To just do trim and alignment to genome write -s "tr-ge"
	-a	adapter sequence for trim (found automaticly if not given), also you can write -a "disable",
		to disable it
	-t	trim front (default 0) How many bases to pre trim reads 5' end,
	        as it frequently represents an untemplated addition during reverse transcription.
	-z	trim tail (default 0) How many bases to pre trim reads on 3' end.
	-A	Alignment type: (default Local, EndToEnd (Local is Local, EndToEnd is force Global))
  -B Discover novel junctions (1, default; 0 keeps only indexed junctions)
  -b Enable fastp base correction for overlapping paired reads
  -M  Max multimapping (default 10) Set to 1 to get only unique reads. Only applies for genome
      step, not the depletion step.
	Path arguments:
	-S   path to STAR (default: ~/bin/STAR-2.7.4a/bin/Linux_x86_64/STAR)
	-P   path to fastp (trimmer) (default: ~/bin/fastp)
	-C	 path to cleaning script, internal

	Less important options:
	-r	resume?: a character (defualt n) (n for new start fresh with file f from point s,
			             (if you want a continue from crash specify the step you want to start
				      from, as defined in -s, start on genome, do -r "ge")
	-m	max cpus allowed (defualt 90)
	-i	include subfolders (defualt n, for no), if you want subfolder do y, for yes.
	-q	a character, Do quality filtering:
	    yes: "default" no: "disable". Uses default fastp QF.
	-k	a character, Keep loaded genomes STAR index:
	    yes (y), no: Remove loaded (n), no shared genome (noShared),
	    default (n)
	-K  Keep contaminant aligned files. Default: "no", alternative: "yes".
	-X  Kept contaminant output file type ("bam", "fastq")
	-u  Keep unaligned reads from genome alignment step.
	      Default "None" (no), alternative: "Fastx""
	-h	this help message

fastp location must be: ~/bin/fastp
STAR location must be: ~/bin/STAR-2.7.4a/bin/Linux_x86_64/STAR

NOTE: if STAR is stuck on load, run this line:
STAR --genomeDir /export/valenfs/data/references/Zv10_zebrafish_allsteps --genomeLoad Remove

example usage: RNA_Align_pipeline_folder.sh -f <in.fastq.gz> -o <out_dir>

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
allow_introns=1
adapter="auto"
quality_filtering="disable"
maxCPU=90
multimap=10
subfolders="n"
trim_front=0
trim_tail=0
paired="no"
align_single=""
cleaning=""
keepLast="n"
keepContam="no"
keepContamType="bam"
keep_unmapped_genome="None"
verbose=1
STAR="$HOME/bin/STAR-2.7.4a/bin/Linux_x86_64/STAR"
fastp="$HOME/bin/fastp"
base_correction=0
in_file_two=""
while getopts ":bvf:o:p:l:T:g:s:a:t:A:B:r:m:k:K:M:S:i:P:I:X:C:q:u:z:h" opt; do
    case $opt in
    b) base_correction=1 ;;
    v)
      verbose=0
      ;;
    f)
        in_dir=$OPTARG
	      ;;
    o)
        out_dir=$OPTARG
        ;;
    p)
	      paired=$OPTARG
        ;;
    l)
        min_length=$OPTARG
        ;;
    T)
        mismatches=$OPTARG
        ;;
    g)
        gen_dir=$OPTARG
        ;;
    s)
        steps=$OPTARG
        ;;
    a)
        adapter=$OPTARG
        ;;
    q)
	      quality_filtering=$OPTARG
        ;;
    t)
        trim_front=$OPTARG
        ;;
    z)
        trim_tail=$OPTARG
        ;;
    A)
        alignment=$OPTARG
        ;;
    B)
        allow_introns=$OPTARG
        ;;
    r)
      	resume=$OPTARG
        ;;
    m)
      	maxCPU=$OPTARG
        ;;
    M)
      	multimap=$OPTARG
        ;;
    i)
      	subfolders=$OPTARG
        ;;
    S)
      	STAR=$OPTARG
        ;;
    P)
      	fastp=$OPTARG
        ;;
    C)
      	cleaning=$OPTARG
        ;;
    I)
      	align_single=$OPTARG
        ;;
    K)
      	keepContam=$OPTARG
        ;;
    k)
      	keepLast=$OPTARG
        ;;
    X)
      	keepContamType=$OPTARG
        ;;
    u)
      	keep_unmapped_genome=$OPTARG
        ;;
    h)
        usage
        exit
        ;;
    :|?)
        echo "Invalid option or missing value: -$OPTARG"
        usage
        exit 1
        ;;
    esac
done

fail() { echo "ERROR: $*" >&2; exit 1; }
log() { if (( verbose )); then echo "$@"; fi; }
[[ -d "$in_dir" ]] || fail "Input directory does not exist: $in_dir"
[[ -n "$out_dir" ]] || fail "Output directory (-o) is required."
[[ -f "$align_single" ]] || fail "Single-library script not found: $align_single"
case "$paired" in yes|no) ;; *) fail "Paired end mode must be yes or no." ;; esac
if [[ "$steps" == all ]]; then
  if [[ -d "$gen_dir/contaminants_genomeDir" ]]; then steps=tr-co-ge
  else steps=tr-ph-rR-nc-tR-ge; fi
fi
IFS='-' read -r -a steps_array <<< "$steps"
[[ "$resume" == n || "-$steps-" == *"-$resume-"* ]] || fail "Resume step is not in steps: $resume"
# Find regular files and symlinks, keeping paths intact. Sorted adjacent files form pairs.
find_args=("$in_dir")
[[ "$subfolders" != n ]] || find_args+=(-maxdepth 1)
files=()
while IFS= read -r -d '' file; do
  [[ "$file" =~ \.(fasta|fa|fastq|fq)(\.gz)?$ ]] && files+=("$file")
done < <(find -L "${find_args[@]}" -type f -print0 | sort -z)
count=${#files[@]}
(( count > 0 )) || fail "No FASTA/FASTQ files found in $in_dir"
stride=1
if [[ "$paired" == yes ]]; then
  (( count % 2 == 0 )) || fail "Paired end input must contain an even number of files."
  stride=2
fi
# The output names use basenames: reject collisions before processing any library.
basenames=()
for ((i=0; i<count; i+=stride)); do
  name=$(basename "${files[i]}"); name=${name%.gz}; name=${name%.*}
  for seen in "${basenames[@]}"; do
    [[ "$name" != "$seen" ]] || fail "Duplicate library basename: $name"
  done
  basenames+=("$name")
done
mkdir -p "$out_dir"
log "Processing $((count / stride)) libraries; steps: $steps"
common=(-o "$out_dir" -a "$adapter" -q "$quality_filtering" -s "$steps"
        -l "$min_length" -T "$mismatches" -g "$gen_dir" -m "$maxCPU" -M "$multimap"
        -A "$alignment" -B "$allow_introns" -t "$trim_front" -z "$trim_tail"
        -K "$keepContam" -X "$keepContamType" -u "$keep_unmapped_genome" -P "$fastp" -S "$STAR")
(( ! base_correction )) || common+=(-b)
(( verbose )) || common+=(-v)
started=0
[[ "$resume" != n ]] || started=1
for current in "${steps_array[@]}"; do
  [[ "$current" != "$resume" ]] || started=1
  (( started )) || continue
  for ((i=0; i<count; i+=stride)); do
    keep=y
    # noShared must apply to every library, not just the last library.
    if [[ "$keepLast" == noShared ]] || (( i + stride == count )); then keep=$keepLast; fi
    log "Step $current; library $((i / stride + 1))/$((count / stride))"
    cmd=(bash "$align_single" "${common[@]}" -f "${files[i]}" -r "$current" -k "$keep")
    [[ "$paired" != yes ]] || cmd+=(-F "${files[i+1]}")
    "${cmd[@]}"
  done
done
if [[ -n "$cleaning" ]]; then bash "$cleaning" "$out_dir"; fi
printf '%q ' "$0" "$@" > "$out_dir/runCommand.log"
printf '\n' >> "$out_dir/runCommand.log"
log "All requested steps completed."
