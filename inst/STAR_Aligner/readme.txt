ORFik STAR / fastp pipeline
==========================

Use STAR.index(), STAR.align.single(), or STAR.align.folder() from R. The
bundled Bash scripts implement these functions; run a script with -h to see
its command-line options. Linux, macOS and WSL are supported.

Examples (paths are placeholders):

  STAR.align.folder("reads", "processed", "references/STAR_index",
                    steps = "tr-co-ge", paired.end = TRUE,
                    base.correction = TRUE)

  STAR.align.single("reads/sample.fastq.gz", output.dir = "processed",
                    index.dir = "references/STAR_index", steps = "tr-ge",
                    allow.introns = FALSE)

Splice junctions
---------------

allow.introns = TRUE (default) passes --alignIntronMax 0 to STAR: discover
novel junctions, using STAR's automatic maximum intron length.
allow.introns = FALSE passes --alignIntronMax 1: suppress novel junctions.
Junctions supplied in the STAR index can still align in either mode.

Earlier affected ORFik versions passed 1 even when TRUE was requested.
The corrected default can change alignment results. Set FALSE explicitly
when you want to preserve the previous indexed-junction-only behavior.

SJ.out.tab contains both annotated and unannotated junction coordinates.
Column 6 is 1 for junctions annotated in STAR's database and 0 otherwise;
the existence of this file alone does not indicate novel junction discovery.

Paired reads and base correction
-------------------------------

Folder mode pairs adjacent files in sorted filename order. Ensure the
mates sort together, for example sample_1.fastq.gz and sample_2.fastq.gz.
Mate identity is not inferred from FASTQ records. Mixed single/paired
libraries must be processed separately. Duplicate output basenames in
recursive folder input are rejected.

base.correction = TRUE enables fastp --correction during trimming. It
corrects mismatches in overlapping paired reads using base quality; it
does not merge pairs. The default is FALSE. It requires paired input and
tr in steps. It has no effect when resuming after trimming.

Steps and logs
--------------

Specify steps in processing order, separated by hyphens:
tr (fastp), co (merged contaminants), ph (PhiX), rR (rRNA), nc (ncRNA),
tR (tRNA), ge (genome). Choose co or individual depletion steps, not both.
"all" selects tr-co-ge if a merged contaminant index exists, otherwise
tr-ph-rR-nc-tR-ge. Indices must exist for the selected alignment steps.

Each step reports the files it actually reads, the output paths, and the
exact tool command. verbose = FALSE hides ORFik's progress messages;
STAR and fastp diagnostics remain visible. Tool failures stop processing
and produce an R error. Folder mode collects STAR logs and junction tables
under each stage's LOGS directory after successful completion.

For resume, keep the original steps and specify the step to resume from.
STAR.align.single runs only that step; STAR.align.folder runs it and all
subsequent steps. Existing intermediate files must be present.

keep.index.in.memory = "noShared" disables shared-memory indices for every
library. TRUE keeps indices loaded; FALSE removes them after the final
library in each folder stage.

References
----------
STAR manual: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
Annotated junction exception: https://github.com/alexdobin/STAR/issues/180
fastp correction: https://github.com/OpenGene/fastp#base-correction-for-pe-data
