#!/bin/bash
# Stop immediately on failed processing commands, including pipeline failures.
set -eo pipefail

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
		Write your wanted steps, seperated by "-". Use processing order, with trimming first and genome alignment last.
		To just do trim and alignment to genome write -s "tr-ge"
	-a	adapter sequence for trim (found automaticly if not given), also you can write -a "disable",
		to disable it or "standard" to get "AAAAAAAAAA", the illumina standard sequence.
	-t	trim front (default 0) How many bases to pre trim reads on 5' end,
	        as it frequently represents an untemplated addition during reverse transcription.
	-z	trim tail (default 0) How many bases to pre trim reads on 3' end.
	-A	Alignment type: (default Local, EndToEnd (Local is Local, EndToEnd is force Global))
  -B Discover novel junctions (1, default; 0 keeps only indexed junctions)
  -b Enable fastp base correction for overlapping paired reads
	Path arguments:
	-S      path to STAR (default: ~/bin/STAR-2.7.4a/bin/Linux_x86_64/STAR)
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
	-K  Keep contaminant aligned files. Default: "no", alternative: "yes".
	-X  Kept contaminant output file type ("bam", "fastq")
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

# Pipeline progress is optional; tool diagnostics remain visible.
log() {
  if (( verbose == 1 )); then echo "$@"; fi
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
trim_front=0
trim_tail=0
keep="n"
keepContam="no"
keepContamType="bam"
keep_unmapped_genome="None"
verbose=1
STAR="$HOME/bin/STAR-2.7.4a/bin/Linux_x86_64/STAR"
fastp="$HOME/bin/fastp"
base_correction=0
in_file_two=""
while getopts ":bvf:F:o:l:T:g:s:a:t:A:B:r:m:M:K:k:p:S:P:X:q:u:z:h" opt; do
    case $opt in
    b) base_correction=1 ;;
    v)
      verbose=0
      ;;
    f)
        in_file=$OPTARG
	      ;;
    F)
      	in_file_two=$OPTARG
	      ;;
    o)
        out_dir=$OPTARG
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
    S)
      	STAR=$OPTARG
        ;;
    P)
      	fastp=$OPTARG
        ;;
    k)
      	keep=$OPTARG
        ;;
    K)
      	keepContam=$OPTARG
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
[[ -n "$out_dir" ]] || fail "Output directory (-o) is required."
[[ -f "$in_file" ]] || fail "Input file does not exist: $in_file"
[[ -z "$in_file_two" || -f "$in_file_two" ]] || fail "Read 2 does not exist: $in_file_two"
case "$allow_introns" in
  1|TRUE|true) intron_max=0 ;;
  0|FALSE|false) intron_max=1 ;;
  *) fail "Allow introns (-B) must be 1/TRUE or 0/FALSE." ;;
esac
case "$keep" in
  y) genome_load=LoadAndKeep ;;
  n) genome_load=LoadAndRemove ;;
  noShared) genome_load=NoSharedMemory ;;
  *) fail "Index memory mode (-k) must be y, n or noShared." ;;
esac
[[ "$keepContamType" == bam ]] || fail "Only BAM contaminant output is supported."
if [[ "$steps" == all ]]; then
  if [[ -d "$gen_dir/contaminants_genomeDir" ]]; then steps=tr-co-ge
  else steps=tr-ph-rR-nc-tR-ge; fi
fi
IFS='-' read -r -a steps_array <<< "$steps"
[[ -n "$steps" && "$steps" != *- ]] || fail "Empty processing step."
previous_rank=0
for step in "${steps_array[@]}"; do
  case "$step" in
    tr) rank=1 ;; co) rank=2 ;; ph) rank=3 ;; rR) rank=4 ;;
    nc) rank=5 ;; tR) rank=6 ;; ge) rank=7 ;;
    *) fail "Unknown processing step: $step" ;;
  esac
  (( rank > previous_rank )) || fail "Steps must be unique and in processing order: tr-co-ph-rR-nc-tR-ge."
  previous_rank=$rank
  if [[ "-$steps-" == *-co-* && "$step" =~ ^(ph|rR|nc|tR)$ ]]; then
    fail "Use merged contaminants (co) or individual depletion steps, not both."
  fi
done
[[ "$resume" == n || "-$steps-" == *"-$resume-"* ]] || fail "Resume step is not in steps: $resume"
if (( base_correction )); then
  [[ -n "$in_file_two" ]] || fail "Base correction requires paired end reads."
  [[ "-$steps-" == *-tr-* ]] || fail "Base correction requires trimming (tr)."
fi
mkdir -p "$out_dir"
ibn=$(basename "$in_file")
ibn=${ibn%.gz}; ibn=${ibn%.*}

# Populate an array, preserving spaces and shell metacharacters in file paths.
step_outputs() {
  case "$1" in
    tr) outputs=("$out_dir/trim/trimmed_${ibn}.fastq")
        [[ -z "$in_file_two" ]] || outputs+=("$out_dir/trim/trimmed2_${ibn}.fastq") ;;
    co) prefix="$out_dir/contaminants_depletion/contaminants_${ibn}_" ;;
    ph) prefix="$out_dir/phix_depletion/PhiX_${ibn}_" ;;
    rR) prefix="$out_dir/rRNA_depletion/rRNA_${ibn}_" ;;
    nc) prefix="$out_dir/ncRNA_depletion/ncRNA_${ibn}_" ;;
    tR) prefix="$out_dir/tRNA_depletion/tRNA_${ibn}_" ;;
    ge) prefix="$out_dir/aligned/${ibn}_" ;;
  esac
  if [[ "$1" != tr ]]; then
    outputs=("${prefix}Unmapped.out.mate1")
    [[ -z "$in_file_two" ]] || outputs+=("${prefix}Unmapped.out.mate2")
  fi
  return 0
}
run_command() {
  if (( verbose )); then
    printf '  Command:'; printf ' %q' "$@"; printf '\n'
  fi
  local status=0
  "$@" || status=$?
  if (( status != 0 )); then
    echo "ERROR: $label failed for $ibn (exit $status). See the tool output above." >&2
    exit "$status"
  fi
}
inputs=("$in_file")
[[ -z "$in_file_two" ]] || inputs+=("$in_file_two")
for step in "${steps_array[@]}"; do
  step_outputs "$step"
  if [[ "$resume" == n || "$resume" == "$step" ]]; then
    for input in "${inputs[@]}"; do
      [[ -f "$input" ]] || fail "Missing input for $step: $input"
    done
    case "$step" in
      tr) label="Trimming (fastp)" ;;
      co) label="Contaminant depletion"; index=contaminants_genomeDir ;;
      ph) label="PhiX depletion"; index=PhiX_genomeDir ;;
      rR) label="rRNA depletion"; index=rRNA_genomeDir ;;
      nc) label="ncRNA depletion"; index=ncRNA_genomeDir ;;
      tR) label="tRNA depletion"; index=tRNA_genomeDir ;;
      ge) label="Genome alignment"; index=genomeDir ;;
    esac
    log "[$step] $label: $ibn"
    if (( verbose )); then printf '  Input: %s\n' "${inputs[@]}"; fi
    if [[ "$step" == tr ]]; then
      mkdir -p "$out_dir/trim"
      if (( verbose )); then printf '  Output: %s\n' "${outputs[@]}"; fi
      cmd=("$fastp" --in1 "${inputs[0]}" --out1 "${outputs[0]}"
           --json "$out_dir/trim/report_${ibn}.json" --html "$out_dir/trim/report_${ibn}.html"
           --trim_front1 "$trim_front" --trim_tail1 "$trim_tail"
           --length_required "$min_length" --thread "$((maxCPU < 16 ? maxCPU : 16))")
      if [[ -n "$in_file_two" ]]; then
        cmd+=(--in2 "${inputs[1]}" --out2 "${outputs[1]}"
              --trim_front2 "$trim_front" --trim_tail2 "$trim_tail")
      fi
      (( ! base_correction )) || cmd+=(--correction)
      [[ "$quality_filtering" != disable ]] || cmd+=(--disable_quality_filtering)
      case "$adapter" in
        auto|"") ;;
        autoPE) cmd+=(--detect_adapter_for_pe) ;;
        disable) cmd+=(--disable_adapter_trimming) ;;
        *) case "$adapter" in
             standard) adapter=AAAAAAAAAA ;; illumina) adapter=AGATCGGAAGAGC ;;
             small_RNA) adapter=TGGAATTCTCGG ;; nextera) adapter=CTGTCTCTTATA ;;
             ingolia12) adapter=CTGTAGGCACCATCAAT ;;
           esac
           cmd+=(--adapter_sequence "$adapter") ;;
      esac
    else
      [[ -d "$gen_dir/$index" ]] || fail "Missing STAR index for $step: $gen_dir/$index"
      mkdir -p "$(dirname "$prefix")"
      reader=-
      [[ "${inputs[0]}" != *.gz ]] || reader=zcat
      if [[ ${#inputs[@]} == 2 ]]; then
        if [[ "$reader" == zcat && "${inputs[1]}" != *.gz || "$reader" == - && "${inputs[1]}" == *.gz ]]; then
          fail "Both mates must use the same compression for STAR."
        fi
      fi
      cmd=("$STAR" --readFilesIn "${inputs[@]}" --genomeDir "$gen_dir/$index"
           --genomeLoad "$genome_load" --outFileNamePrefix "$prefix"
           --outFilterMatchNmin "$min_length" --readFilesCommand "$reader"
           --limitIObufferSize 50000000)
      if [[ "$step" == ge ]]; then
        log "  Output: ${prefix}Aligned.sortedByCoord.out.bam"
        if [[ "$intron_max" == 0 ]]; then
          log "  Splicing: indexed and novel junctions (automatic maximum intron length)"
        else
          log "  Splicing: indexed junctions only (novel junction discovery disabled)"
        fi
        cmd+=(--outSAMtype BAM SortedByCoordinate --outReadsUnmapped "$keep_unmapped_genome"
              --runThreadN "$((maxCPU < 80 ? maxCPU : 80))" --limitBAMsortRAM 30000000000
              --alignEndsType "$alignment" --alignIntronMax "$intron_max"
              --outFilterMultimapNmax "$multimap" --outFilterMismatchNmax "$mismatches")
      else
        if (( verbose )); then printf '  Unmapped output: %s\n' "${outputs[@]}"; fi
        case "$step" in co|rR) thread_cap=90 ;; nc) thread_cap=80 ;; *) thread_cap=70 ;; esac
        cmd+=(--outReadsUnmapped Fastx --runThreadN "$((maxCPU < thread_cap ? maxCPU : thread_cap))")
        if [[ "$keepContam" == yes ]]; then
          cmd+=(--outSAMtype BAM Unsorted --outSAMmode Full)
          log "  Contaminant BAM: ${prefix}Aligned.out.bam"
        else cmd+=(--outSAMtype None --outSAMmode None); fi
        # Preserve the existing depletion settings; these do not use allow.introns.
        case "$step" in ph|rR|tR) cmd+=(--alignIntronMax 1) ;; esac
        [[ "$step" != tR ]] || cmd+=(--seedPerWindowNmax 20 --outFilterMultimapNmax 20)
      fi
    fi
    run_command "${cmd[@]}"
  fi
  inputs=("${outputs[@]}")
done
