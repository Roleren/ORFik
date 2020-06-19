Welcome to STAR aligner and indexing script

There are 3 main scripts you can use:

1. Make star Index, file: STAR_MAKE_INDEX.sh
2. Align data from specified file: RNA_Align_pipeline.sh
3. Align all data from specified folder: RNA_Align_pipeline_folder.sh

Arguments for Star Index creation:
OPTIONS:
	-o	   output folder for all indices
	-s   	   species to use (zebrafish, soon human and yeast)
	-phix 	   path to phix fasta/fasta.gz file
	-rRNA	   path to rrna fasta/fasta.gz file
	-ncRNA	   path to ncRNA fasta/fasta.gz file
	-tRNA	   path to trna fasta/fasta.gz file
	-genome    path to species whole genome fasta/fasta.gz file
	-genomeGTF path to species gtf file
	-outGenome path to use as output for species genome, if other than rest
	-h	   this help message

Arguments for aligning data:
OPTIONS:
	-f	path to input fasta file. Also define -F if paired end! Must be file type of: fasta, fa, fastq, fq or gz
	-F 	path to input fastq file 2 (paired end reads 2) Must be file type of: fasta, fa, fastq, fq or gz
	-o	path to output dir
	-l	minimum length of reads (default: 15)
	-g	genome dir for all indices (Standard is zebrafish: danrerio10, change to human index if needed etc)
	-s	steps of depletion and alignment wanted:
		(a string: which steps to do? (default: "tr-ge", write "all" to get all: "tr-ph-rR-nc-tR-ge")
			 tr: trimming, ph: phix depletion, rR: rrna depletion, 
			 nc: ncrna depletion, tR: trna depletion, ge: genome alignment) 
		Write your wanted steps, seperated by "-". Order does not matter.
		To just do trim and alignment to genome write -s "tr-ge"
	-a	adapter sequence for trim (found automaticly if not given), also you can write -a "disable", 
		to disable it or "standard" to get "AAAAAAAAAA", the illumina standard sequence for 5' end.
	-t	trim front (default 3) How many bases to pre trim 5′ end of each read, 
	        as it frequently represents an untemplated addition during reverse transcription.
	-A	Alignment type: (default Local, EndToEnd (Local is Local, EndToEnd is force Global))

	Less important options:
	-r	resume?: a character (defualt n) (n for new start fresh with file f from point s, 
			             (if you want a continue from crash specify the step you want to start 
				      from, as defined in -s, start on genome, do -r "ge")
	-m	max cpus allowed (defualt 90)
	-h	this help message

NOTE: When running on folder, f is folder directory not a file! F does not exist for folder, as it auto detects pairs of data.
Also, remember you need to run the alignment script from the folder with the fasta files!
Only works on unix systems, not windows.

