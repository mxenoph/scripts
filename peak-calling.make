
# Pseudoreps # {{{

.PHONY: make-pseudoreps
## Creating pseudoreplicates for IDR framework
make-pseudoreps: ${pseudoreps}

${IDR_path}/pseudoreps/%_pseudoreps00.bam ${IDR_path}/pseudoreps/%_pseudoreps01.bam:: ${map_path}/%.bam
	${bsub} -M ${memlimit} -R "rusage[mem=2000]" make-pseudoreps.sh -b $< -o $(@D)

${IDR_path}/pseudoreps/%_pseudoreps00.bam ${IDR_path}/pseudoreps/%_pseudoreps01.bam:: ${map_path}/downsample/%.bam
	${bsub} -M ${memlimit} -R "rusage[mem=2000]" make-pseudoreps.sh -b $< -o $(@D)
# }}}

# Peak Calling # {{{
.PHONY: peak-calling
## Calling broad and narrow peaks for pooled data and independent replicates adn relaxed peaks for the IDR framework
peak-calling: pool-reps ${peaks_broad} ${peaks_sharp} make-pseudoreps ${peaks_relaxed_IDR} ${peaks_pseudoreps}

# For each replicate call peaks independently. Do not split target_ctrl into 2 files, it returns empty
${macs_path_broad}%_peaks.broadPeak: ${map_path}%.bam
	$(eval target_ctrl := $(call format-ctrl,$(call get-ctrl,${map_path},$(call get-base,$*),$(call get-ctrl-type,$*),${input})))
	${bsub} -n 8 -M ${memlimit} -R 'rusage[mem=9000]'\
		"run-macs2.sh -i $< $(target_ctrl) -g mm -o $(patsubst %/, %, $(dir $(@D))) -a ${assembly} -b"

${macs_path_sharp}%_peaks.narrowPeak: ${map_path}%.bam
	$(eval target_ctrl := $(call format-ctrl,$(call get-ctrl,${map_path},$(call get-base,$*),$(call get-ctrl-type,$*),${input})))
	${bsub} -n 16 -M ${memlimit} -R 'rusage[mem=9000]'\
		"run-macs2.sh -i $< $(target_ctrl) -g mm -o $(patsubst %/, %, $(dir $(@D))) -a ${assembly}"# }}}

# Relaxed peaks for IDR # {{{

${IDR_path}/macs2/sharp%_peaks.narrowPeak:: ${map_path}%.bam
	$(eval target_ctrl := $(call get-ctrl,${map_path},$(call get-base,$*),$(call get-ctrl-type,$*),${input}))
	${bsub} -n 16 -M ${memlimit} -R 'rusage[mem=9000]'\
		"run-macs2.sh -i $< -c $(target_ctrl) -g mm -o $(patsubst %/, %, $(dir $(@D))) -a ${assembly} -r"

${IDR_path}/macs2/sharp%_peaks.narrowPeak:: ${map_path}/downsample%.bam
	$(eval target_ctrl := $(call get-ctrl,${map_path},$(call get-base,$*),$(call get-ctrl-type,$*),${input}))
	${bsub} -n 16 -M ${memlimit} -R 'rusage[mem=9000]'\
		"run-macs2.sh -i $< -c $(target_ctrl) -g mm -o $(patsubst %/, %, $(dir $(@D))) -a ${assembly} -r"

${IDR_path}/macs2/sharp%_peaks.narrowPeak:: ${IDR_path}/pseudoreps%.bam
	$(eval target_ctrl := $(call get-ctrl,${map_path},$(call get-base,$*),$(call get-ctrl-type,$*),${input}))
	${bsub} -n 16 -M ${memlimit} -R 'rusage[mem=9000]'\
		"run-macs2.sh -i $< -c $(target_ctrl) -g mm -o $(patsubst %/, %, $(dir $(@D))) -a ${assembly} -r"

# }}}

# Conservative peaks from IDR framework# {{{

.PHONY: get-IDR-threshold
get-IDR-threshold: ${IDR_path}/IDR-thresholds.tsv ${IDR_path}/self-consistency-IDR-thresholds.tsv ${IDR_path}/pooled-consistency-IDR-thresholds.tsv

${IDR_path}/IDR-thresholds.tsv: $(call FILTER_OUT,pooled,${peaks_relaxed_IDR})
	wc -l $^ | sed 's/ \+/\t/g' > ${IDR_path}/number-pre-IDR-peaks.txt
	${bsub} "Rscript ~/source/get-IDR-threshold.R -p ${IDR_path}/number-pre-IDR-peaks.txt -o $@"

${IDR_path}/self-consistency-IDR-thresholds.tsv: $(call FILTER_OUT,pooled,${peaks_pseudoreps})
	wc -l $^ | sed 's/ \+/\t/g' > ${IDR_path}/ps-number-pre-IDR-peaks.txt
	${bsub} "Rscript ~/source/get-IDR-threshold.R -p ${IDR_path}/ps-number-pre-IDR-peaks.txt -o $@"

${IDR_path}/pooled-consistency-IDR-thresholds.tsv: $(call KEEP,pooled,${peaks_pseudoreps})
	wc -l $^ | sed 's/ \+/\t/g' > ${IDR_path}/pooled-ps-number-pre-IDR-peaks.txt
	${bsub} "Rscript ~/source/get-IDR-threshold.R -p ${IDR_path}/pooled-ps-number-pre-IDR-peaks.txt -o $@ -t"
# }}}

# Can't remember why I used secondary expansion here initially but messes up the macros defined below if used here
#SECONDEXPANSION:

# Macros and function specific to the IDR analysis# {{{

get-relaxed = $(call FILTER_OUT,pooled,$(call KEEP,$(1),${peaks_relaxed_IDR}))
get-relaxed-ps = $(call FILTER_OUT,pooled,$(call KEEP,$(1),${peaks_pseudoreps}))
get-relaxed-pool = $(call KEEP,pooled,$(call KEEP,$(1),${peaks_pseudoreps}))
get-exp-thr = $(shell grep $(1) $(2) | cut -f 1 | sed 's/.*_//' | sed 's/.narrowPeak//')
get-top = $(shell grep $(1) $(2) | cut -f 4)

define run-IDR
${IDR_path}/final-set/$(1)_conservative_$(5).narrowPeak: ${IDR_path}/macs2/sharp/$(1)_pooled_peaks.narrowPeak .FORCE
	${bsub} -n 6 -M ${memlimit} -R "rusage[mem=10000]" "run-IDR.sh -p '$(2)' -o ${IDR_path} -n '0 F p.value' -d ${IDR_path}/macs2/sharp/$(1)_pooled_peaks.narrowPeak -t $(6) -q $(5) -e '$(3)' -s '$(4)' -a $(7) -b $(9) -f $(8) -g $(10)"
	#run-IDR.sh -p '$(2)' -o ${IDR_path} -n '0 F p.value' -d ${IDR_path}/macs2/sharp/$(1)_pooled_peaks.narrowPeak -t $(6) -q $(5) -e '$(3)' -s '$(4)' -a $(7) -b $(9) -f $(8) -g $(10)
endef

# }}}

# Implementations of the IDR method by ENCODE
.PHONY: IDR
## IDR framework
IDR: ${peaks_relaxed_IDR} ${peaks_pseudoreps} get-IDR-threshold ${IDR_path}/plots/flag_replicates.pdf

IDR_conservative = $(word 5,$(addprefix ${IDR_path}/final-set/, $(shell tail -n +2 ${IDR_path}/IDR-thresholds.tsv | cut -f 1 )))

# To make it readable this is split with " \" (space backslash). Foreach commands need to be split like splitting recipe lines
# (backslash is not replaced by space, thus space needs to be added before) otherwise the arguments after the first backslash will not be passed to run-IDR
$(foreach r,$(notdir ${IDR_reps}), \
	$(eval $(call run-IDR,$(r),$(call get-relaxed,$(r)),$(call get-relaxed-ps,$(r)), \
	$(call get-relaxed-pool,$(r)),$(call get-exp-thr,$(r),${IDR_path}/IDR-thresholds.tsv),$(call get-top,$(r),${IDR_path}/IDR-thresholds.tsv), \
	$(call get-top,$(r),${IDR_path}/self-consistency-IDR-thresholds.tsv),$(call get-exp-thr,$(r),${IDR_path}/self-consistency-IDR-thresholds.tsv), \
	$(call get-top,$(r),${IDR_path}/pooled-consistency-IDR-thresholds.tsv),$(call get-exp-thr,$(r),${IDR_path}/pooled-consistency-IDR-thresholds.tsv))))

${IDR_path}/replicate_thresholds.txt: ${IDR_conservative}
	for f in `ls ${IDR_path}/batch-consistency/*-overlapped-peaks.txt`; do echo $$f; awk -v idr=${IDR_threshold} '$$11 <= idr {print $$0}' $$f | wc -l; done | paste - - > $@

${IDR_path}/plots/flag_replicates.pdf: ${IDR_path}/replicate_thresholds.txt
	Rscript -e "source('~/source/Rscripts/plotting-functions.R'); pdf('$@'); flag_replicates('$<'); dev.off()"

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
