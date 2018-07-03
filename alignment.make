.PHONY: alignment# {{{
## Aligning with bowtie
alignment: ${mapped_reads} #${ddup_reads}

# https://www.gnu.org/software/make/manual/html_node/Double_002dColon.html
${map_path}%.bam:: ${data_path}%.fastq.gz ${data_path}/fastqc/%_fastqc.zip
	${bsub} -n 4 "run-bowtie.sh -f $< -o ${map_path}/ -c ${data_path}/fastqc"
	
# }}}

.PHONY: ddup# {{{
## Removing duplicate reads
ddup: ${ddup_reads}

${map_path}/ddup/%_ddup.bam: ${map_path}/%.bam
	${bsub} -n 4 "${java} ${picard}/MarkDuplicates.jar INPUT=$< OUTPUT=$@ ASSUME_SORTED=true REMOVE_DUPLICATES=true M=$(strip $(patsubst %.bam, %.log, $@))"

# }}}

.PHONY: downsample# {{{
## Downsampling to match sequencing depth in both bam and ddup.bam
downsample: ${subset} ${ddup_subset}

${map_path}/downsample/%_downsample.bam: ${map_path}/%.bam ${map_path}/alignment-summarised-stats.tsv
	$(eval perc := $(shell grep $(patsubst %_downsample.bam, %, $(@F)) ${map_path}/alignment-summarised-stats.tsv | cut -f 12 | sed 's/.*\.//'))
	${bsub} -n 4 "samtools view -h -s 5.$(perc) $< -bo $@"

${map_path}/downsample/ddup/%_downsample.bam: ${map_path}/ddup/%.bam ${map_path}/alignment-summarised-stats.tsv
	$(eval perc := $(shell grep $(patsubst %_ddup_downsample.bam, %, $(@F)) ${map_path}/alignment-summarised-stats.tsv | cut -f 9 | sed 's/.*\.//'))
	${bsub} -n 4 "samtools view -h -s 5.$(perc) $< -bo $@"# }}}

.PHONY:sam-to-bam# {{{
## Converting sam to bam
sam-to-bam: ${mapped_reads}

${map_path}%.bam:: ${map_path}%.sam
	${bsub} -M 8000 -n 8 -R 'rusage[mem=8000]'  "samtools view -bS $< | samtools sort - $(patsubst %.bam, %, $@)"# }}}

.PHONY:index# {{{
## Indexing bam files (including deduplicated)
index: ${indexed_reads} ${ddup_indexed}

# Paired end indexing throws an error
${map_path}%.bam.bai: ${map_path}%.bam
	${bsub} samtools index $<

${map_path}/ddup%.bam.bai: ${map_path}/ddup%.bam
	${bsub} samtools index $<
	
# }}}

.PHONY: mapping-stats# {{{
## Alignment statistics
mapping-stats: ${mapping_stats} ${map_path}/bowtie-alignment-statistics.tsv

${map_path}%.stats: ${map_path}%.bam
	${bsub} "samtools stats $< > $@"

${map_path}%.n_mapped: ${map_path}%.bam
	${bsub} "samtools idxstats $< > $@"

${map_path}/bowtie-alignment-statistics.tsv: ${mapped_reads}
	${bsub} "~/source/count-bowtie-reads.sh -b $(<D)"
	${bsub} "Rscript ~/source/plot-alignment-stats.R -f $@ -o $(@D)"
	
# }}}

.PHONY: coverage# {{{
## Calculate coverage on bam and ddup.bam (simple scaling on total mapped reads)
coverage: pool-reps ${bigwig}

${coverage_path}%.bw: ${map_path}%.bam
	${bsub} -M ${memlimit} "genomeCoverageBed -bg -ibam $< -split -g ${chromosomes} > $(@D)$*.bedgraph"
	${bsub} -M ${memlimit} "sort -k1,1 -k2,2n $(@D)$*.bedgraph > $(@D)$*.sorted.bedgraph"
	mv $(@D)$*.sorted.bedgraph $(@D)$*.bedgraph
	${bsub} -M ${memlimit} "bedGraphToBigWig $(@D)$*.bedgraph ${chromosomes} $@"
	rm $(@D)$*.bedgraph

${coverage_path}/ddup%.bw: ${map_path}/ddup%.bam
	${bsub} -M ${memlimit} "genomeCoverageBed -bg -ibam $< -split -g ${chromosomes} > $(@D)$*.bedgraph"
	${bsub} -M ${memlimit} "sort -k1,1 -k2,2n $(@D)$*.bedgraph > $(@D)$*.sorted.bedgraph"
	mv $(@D)$*.sorted.bedgraph $(@D)$*.bedgraph
	${bsub} -M ${memlimit} "bedGraphToBigWig $(@D)$*.bedgraph ${chromosomes} $@"
	rm $(@D)$*.bedgraph

# }}}

#.PHONY: bowtie2# {{{
## Running bowtie2 instead of bowtie for ChIP-RX paper protocol
#bowtie2: ${mapped_reads}

# https://www.gnu.org/software/make/manual/html_node/Double_002dColon.html
#${map_path}%.bam:: ${data_path}%.fastq.gz ${data_path}/fastqc/%_fastqc.zip
#	echo ${bsub} -n 4 "run-bowtie2.sh -f $< -g ${genome} -o ${map_path}/ -c ${data_path}/fastqc"# }}}


# vim: ft=make
