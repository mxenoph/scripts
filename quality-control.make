# FASTQ QC# {{{
.PHONY: fastq-qc
fastq-qc: ${fastqc_zip}

${data_path}/fastqc/%_fastqc.zip: ${data_path}/%.fastq.gz
	${bsub} -n 6 -R 'rusage[mem=6000]' "run-fastqc.sh $<"

# }}}

# ChIP quality control# {{{
.PHONY: chip-qc
## Calculates cross-correlation and makes custom fingerprint plots to investigate enrichment
chip-qc: ${cross_correlation} ${enrichment_files} ${enrichment_plots}

# use 6 nodes, savp creates the cross correlation plot
# run_spp_parallel.R is a modified version of phantompeaktools/run_spp.R to use 
# makeForkCluster from parallel package (instead of makeCluster from snow, which doesn't play well
# with LSF in this case). Similar problem as here: https://stat.ethz.ch/pipermail/r-sig-hpc/2012-December/001568.html

${qc_path}%.cross-correlation:: ${map_path}%.bam
	${bsub} -n 8 -M ${memlimit} -R "select[ncores=8]" "Rscript ~/local/phantompeakqualtools/run_spp_parallel.R -p=8 -c=$< -savp -tmpdir=${qc_path}/tmp -out=$@"
	mv ${map_path}$*.pdf ${qc_path}

${qc_path}/%.cross-correlation:: ${map_path}/ddup/%.bam
	${bsub} -n 8 -M ${memlimit} -R "select[ncores=8]" "Rscript ~/local/phantompeakqualtools/run_spp_nodups.R -p=8 -c=$< -savp -tmpdir=${qc_path}/tmp -out=$@"
	mv ${map_path}/ddup/$*.pdf ${qc_path}

${qc_path}%.enrichment: ${map_path}%.bam
	${bsub} -n 4 "perl ~/source/check-chip-enrichment.pl ${reference}/${assembly}.fa.fai 1000 $< > $@"

${qc_path}/plots/%.png: ${qc_path}/%.enrichment
	$(eval base := $(shell echo $(notdir $*) | sed 's/-.*//g'))
	$(eval index := $(shell echo $(notdir $*) | sed 's/.*_//g'))
	$(eval target_ctrl := $(call format-ctrl,$(call get-ctrl,${map_path},$(call get-base,$*),$(call get-ctrl-type,$*),${input})))
	${bsub} -M ${memlimit} "Rscript ~/source/plot-chip-enrichment.R -b $< $(target_ctrl) -o $(dir $(@D))"

# }}}

# Correlation # {{{
correlation := $(call FILTER_OUT,7G9_EpiSC,\
	$(call FILTER_OUT,Remco,$(addprefix ${coverage_path}/ddup/,$(addsuffix .sample-correlation, ${conditions}))))
correlation_plots := $(patsubst %.sample-correlation, %.pdf, ${correlation})
on_sets := $(patsubst %.sample-correlation, %, ${correlation})

.PHONY: correlation
## Calculate correlation of samples genome wide (not on specific regions)
correlation: ${correlation} ${correlation_plots}

define COR
$(1).sample-correlation: ${bigwig}
	$(eval per_filter := $(call KEEP,ddup,$(call KEEP,$(notdir $(1)),${bigwig})))
	@echo $(per_filter) | tr ' ' '\n' > $(1).tmp
	Rscript -e "a=read.delim('$(1).tmp', stringsAsFactors=F); a=a[[1]]; b=as.data.frame(t(combn(a,2))); write.table(b, '$(1).pairwise-combn', col.names=F, row.names=F, sep=' ', quote=F)"
	rm $(1).sample-correlation
	${bsub} run-bigWigCorrelate-on-pairs.sh $(1).pairwise-combn $(1).sample-correlation
	rm $(1).tmp

$(1).pdf: $(1).sample-correlation
	${bsub} "Rscript ~/source/plot-bw-correlation.R -s $(1).sample-correlation"

endef

$(foreach cn,${on_sets},$(eval $(call COR,$(cn))))

.PHONY: correlate
correlate: ${correlation} ${correlation_plots}

.PHONY: plot-correlation
plot-correlation: ${correlation_plots}

# }}}

# Filtering blacklisted regions# {{{
.PHONY: filter-peaks
## Remove ENCODE blacklisted regions
filter-peaks: ${filtered_peaks}

${macs_path_sharp}/%_filtered_peaks.narrowPeak: ${macs_path_sharp}/%_peaks.narrowPeak ${blacklist}
	${bsub} "bedtools intersect -a $< -b ${blacklist} -v > $@"

${macs_path_broad}/%_filtered_peaks.broadPeak: ${macs_path_broad}/%_peaks.broadPeak ${blacklist}
	${bsub} "bedtools intersect -a $< -b ${blacklist} -v > $@"

${macs_path_broad}/%_filtered_peaks.gappedPeak: ${macs_path_broad}/%_peaks.gappedPeak ${blacklist}
	${bsub} "bedtools intersect -a $< -b ${blacklist} -v > $@"

${IDR_path}/%_filtered_conservative_${IDR_threshold}.narrowPeak: ${IDR_path}/%_conservative_${IDR_threshold}.narrowPeak ${blacklist}
	${bsub} "bedtools intersect -a $< -b ${blacklist} -v > $@"

# }}}
# vim: ft=make