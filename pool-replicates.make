## Function pools replicates (that have not failed) together for each protein and input # {{{
define REPS
$(1)_pooled.bam: #good-experiments.txt
	$(eval mapped = $(strip $(shell cat good-experiments.txt | tr '\n' ' ')))
	$(eval ds := $(shell ls -d ${map_path}/downsample/*))
	$(eval ds := $(call KEEP,$(addsuffix _,$(notdir $(1))),$(ds)))
	$(eval replace := $(patsubst %_downsample.bam, %.bam,$(notdir $(ds))))
	$(eval tmp = $(call KEEP,$(addsuffix _,$(1)),$(mapped)))
	$(foreach d,$(replace),$(eval tmp = $(call FILTER_OUT,$(d),$(tmp))))
	$(eval tmp := $(tmp) $(ds))
	${bsub} -n 4 -R "rusage[mem=8000]" "samtools merge -f $(1)_pooled.bam $(tmp)"

$(1)_pooled_ddup.bam: #good-experiments-ddup.txt
	$(eval dduped = $(strip $(shell cat good-experiments-ddup.txt | tr '\n' ' ')))
	$(eval ds_ddup := $(shell ls -d ${map_path}/downsample/ddup/*))
	$(eval ds_ddup := $(call KEEP,$(addsuffix _,$(notdir $(1))),$(ds_ddup)))
	$(eval replace_ddup := $(patsubst %_downsample.bam, %.bam,$(notdir $(ds_ddup))))
	$(eval tmp_ddup = $(call KEEP,$(addsuffix _,$(1)),$(dduped)))
	$(foreach dd,$(replace_ddup),$(eval tmp_ddup = $(call FILTER_OUT,$(dd),$(tmp_ddup))))
	$(eval tmp_ddup := $(tmp_ddup) $(ds_ddup))
	${bsub} -n 4 -R "rusage[mem=8000]" "samtools merge -f $(1)_pooled_ddup.bam $(tmp_ddup)"
endef
# Prereq list cannot be accessed with $^ in the rules defined in the function 
# and dependencies are not checked. However the existence of good-experiments.txt 
# is checked before hand so this does not cause a problem. TODO: find a more 
# elegant way to deal with pooling
# }}}

# Excluding failed experiments from pooling # {{{
.PHONY: remove-failed
## Exclude failed experiments from pooling based on contents of user managed failed-experiments.tsv
remove-failed: good-experiments.txt good-experiments-ddup.txt

good-experiments.txt: failed-experiments.txt
	$(eval x:= $(strip $(shell cat failed-experiments.txt 2>/dev/null)))
	$(eval tmp = ${mapped_reads})
	$(foreach failed,$(x),$(eval tmp = $(call FILTER_OUT,$(failed),$(tmp))))
	@echo $(tmp) | tr ' ' '\n' > $@

good-experiments-ddup.txt: failed-experiments.txt
	$(eval x:= $(strip $(shell cat failed-experiments.txt 2>/dev/null)))
	$(eval tmp = ${ddup_reads})
	$(foreach failed,$(x),$(eval tmp = $(call FILTER_OUT,$(failed),$(tmp))))
	@echo $(tmp) | tr ' ' '\n' > $@
# }}}

# Call REPS function to create targets
$(foreach r,${reps},$(eval $(call REPS,$(r))))
$(foreach dr,${ddup_reps},$(eval $(call REPS,$(dr))))

.PHONY: pool-reps
## Pooling replicates
# the rule to downsample is required for the Hendrich chip not for all other datasets
#pool-reps: ${mapped_reads} remove-failed downsample ${pooled} ${ddup_pooled}
pool-reps: ${mapped_reads} remove-failed ${pooled} ${ddup_pooled}
# vim: ft=make
