#!/bin/bash
# The commands describe the first steps of the Bulk RNAseq analysis

# Author : Ouardia Ait-Mohamed # 

# tools and hg38 versions #
# GRChg38.v93
# hg39.v93.gtf
# STAR-2.7.3a
# samtools-1.9
# subread-1.5.1

SCRIPT_PATH="PATH_TO/pre_dc_1281_06/Scripts"

STAR_PATH="PATH_TO/STAR"
SAMTOOLS_PATH="PATH_TO/samtools"
FEATURECOUNTS_PATH="PATH_TO/subread"

PATH=$STAR_PATH:$SAMTOOLS_PATH:$FEATURECOUNTS_PATH:$PATH
export PATH 


mkdir -p $TMPDIR

GENOME1="hg38"
GENOME2="plasmid"
OUT_DIR1=PATH_TO/$GENOME1"
OUT_DIR2=PATH_TO/$GENOME2"
in_dir=PATH_TO/Fastq"
genome_dir=PATH_TO/STAR_index"

trap 'rm -f “$ftmp”' EXIT
ftmp=$(mktemp)|| exit 1
find $in_dir -iname *R1*.gz -printf %f'\n' | cut -d'.' -f1 >$ftmp;


cat $ftmp
while read -r line; do
	name=$line
	out_prefix1="$OUT_DIR1/$name"
	out_prefix2="$OUT_DIR2/$name"

	#STAR mapping on human genome
	rm -r $TMPDIR/STAR
STAR --readFilesIn $in_dir/$name.R1.fastq.gz $in_dir/$name.R2.fastq.gz --readFilesCommand zcat --outFileNamePrefix $out_prefix1. --outMultimapperOrder Random --outFilterMismatchNmax 2 --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 --alignIntronMin 20 --alignIntronMax 6000 --genomeDir $genome_dir/$GENOME1 --outSAMtype BAM SortedByCoordinate --outReadsUnmapped Fastx --outTmpDir $TMPDIR/STAR
	# Zip the unmapped fastq files 
gzip $out_prefix1.Unmapped.out.mate1
gzip $out_prefix1.Unmapped.out.mate2
	# Remove chimeric reads
samtools view -hF 0x904 $out_prefix1.Aligned.sortedByCoord.out.bam  | samtools view -Sb > $OUT_DIR1/filtred/$name.bam
        # Create counting files for the human transcripts
featureCounts -s 1 -a hg38.v93.gtf -o $OUT_DIR1/count/$name.txt $OUT_DIR1/filtred/$name.bam	

       #STAR mapping on plamsid genome
STAR --readFilesIn $out_prefix1.Unmapped.out.mate1.gz $out_prefix1.Unmapped.out.mate2.gz --readFilesCommand zcat --outFileNamePrefix $out_prefix2. --outMultimapperOrder Random --outFilterMismatchNmax 2 --alignEndsType EndToEnd --genomeDir $genome_dir/$GENOME2 --outSAMtype BAM SortedByCoordinate --outTmpDir $TMPDIR/STAR
	# Remove chimeric reads
samtools view -hF 0x904 $out_prefix2.Aligned.sortedByCoord.out.bam  | samtools view -Sb > $OUT_DIR2/filtred/$name.bam
done < $ftmp
