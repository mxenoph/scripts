get-fragment-size = $(shell if [[ -z $(strip $(1)) ]] ; then echo ''; else f=$$(awk -F '\t' '{split($$3,fs,",")} END{for(i=1;i<=NF;i++){if(fs[i] > 50) {print fs[i]; break}}}' < $(1)); echo "$${f}"; fi)

# Normalise to input with bamCompare #{{{
.PHONY: bam-compare
bam-compare: index ${count_bigwig} ${ses_bigwig} save-sf


define build-bam-compare-command =
	mkdir -p $(dir $3)
	$(eval target_ctrl := $(call get-ctrl,${map_path}/ddup,$(call get-base,$1),$(call get-ctrl-type,$1),$(addprefix ${map_path}/ddup/,$(patsubst %.bam, %_ddup.bam, $(notdir ${input})))))
	@echo $(target_ctrl)
	# keeping stem (wihtout pooled) to search for the cross correlation files of the replicates making up the pooled data
	$(eval tmp_1 := $(patsubst %_pooled, %, $(notdir $1)))
	$(eval cc := $(call KEEP,ddup,$(call KEEP,$(tmp_1),${cross_correlation})))
	$(eval fragment_size := )
	$(foreach c,$(cc),$(eval fragment_size := $(fragment_size) $(if $(findstring ${output_path},$(c)), $(call get-fragment-size, $(c)))))
	$(eval fragment_size := $(shell echo $(fragment_size) | Rscript -e "f = file('stdin'); a = readLines(f); tmp = unlist(stringr::str_split(a, ' ')); cat(round(median(as.numeric(tmp))))"))
	# need to run as one bash command otherwise make executes in different shells simultaneously and can not get job_id if I cat the file before the job is submitted. 
	# Last if-condition handles empty target files so that the rule is not satisfied unless the target file is non-empty and thus valid
	${bsub} -M ${memlimit} -n 8 "bamCompare -b1 $2 -b2 $(target_ctrl) -o $3 -of bigwig --scaleFactorsMethod $4 --ratio $5 --centerReads -v -p 8 --binSize 10 --extendReads $(fragment_size)" | grep -Eow "[0-9]+" > $(patsubst %.bw, %.run-info.tmp, $3);
endef

${deeptools_path}/read-count/%_RPGC.bw: ${map_path}/ddup/%_ddup.bam ${map_path}/ddup/%_ddup.bam.bai
	$(call build-bam-compare-command,$*,$<,$@,'readCount','subtract' --normalizeTo1x ${effective_genome})

${deeptools_path}/read-count/%_RPKM.bw: ${map_path}/ddup/%_ddup.bam ${map_path}/ddup/%_ddup.bam.bai
	$(call build-bam-compare-command,$*,$<,$@,'readCount','subtract' --normalizeUsingRPKM)

${deeptools_path}/ses-norm/%_RPGC.bw: ${map_path}/ddup/%_ddup.bam ${map_path}/ddup/%_ddup.bam.bai
	$(call build-bam-compare-command,$*,$<,$@,'SES','subtract' --normalizeTo1x ${effective_genome})

${deeptools_path}/ses-norm/%_RPKM.bw: ${map_path}/ddup/%_ddup.bam ${map_path}/ddup/%_ddup.bam.bai
	$(call build-bam-compare-command,$*,$<,$@,'SES','subtract' --normalizeUsingRPKM)
# }}}

# Removing the generation of the size factor files from the build-bam-compare # {{{
# function as it caused problems with the grep one-liner
# running before completion of the bjob. sf files will only be created for those sample that have no size-factor file, have a .run-info.tmp file with a JOBID
# corresponding to a file present in $PWD.logs
curr_sf_files := $(patsubst %.size-factors, %, $(shell find ${deeptools_path} -maxdepth 2 -type f -name "*size-factors"))
bamCompare_JOBID := $(patsubst %.run-info.tmp, %, $(shell find ${deeptools_path} -maxdepth 2 -type f -name "*tmp"))
target_sf_files := $(shell awk 'NR==FNR{a[$$1]++;next;}!($$0 in a)' <(tr ' ' "\n" <<< "${curr_sf_files}") <(tr ' ' "\n" <<< "${bamCompare_JOBID}"))
target_sf_files := $(patsubst %, %.size-factors, $(shell for f in ${target_sf_files}; do log_file=$$PWD/logs/bamCompare-$$(cat $$f.run-info.tmp).log; if [ -e $$log_file ]; then echo $$f; fi; done ))

.PHONY: save-sf
save-sf: ${target_sf_files}
${deeptools_path}/read-count/%.size-factors: ${deeptools_path}/read-count/%.run-info.tmp
	# Removing everything after line /Individual scale factors/
	grep -A 5 'Size factors using ' $$PWD'/logs/bamCompare-'`cat $<`'.log' | sed '/Individual scale factors/q' | sed 's/.*Size/Size/' > $@
	rm $<

${deeptools_path}/ses-norm/%.size-factors: ${deeptools_path}/ses-norm/%.run-info.tmp
	grep -A 5 'Size factors using ' $$PWD'/logs/bamCompare-'`cat $<`'.log' | sed '/Individual scale factors/q' | sed 's/.*Size/Size/' > $@
	rm $<

# }}}

#Functions to compute matrix# {{{
define compute-region-matrix =
	mkdir -p $(sort $(dir $(3)))flat-files; ${bsub} -n 6 "computeMatrix scale-regions -S $(1) -R $(2) --outFileName $(3) --outFileNameMatrix $(sort $(dir $(3)))flat-files/$(notdir $(3)).tsv --outFileSortedRegions $(3)-sorted-by.bed -bs 10 -p 6 $(4)"
endef

define compute-point-matrix =
	mkdir -p $(sort $(dir $(3)))flat-files; ${bsub} -n 6 "computeMatrix reference-point -S $(1) -R $(2) --outFileName $(3) --outFileNameMatrix $(sort $(dir $(3)))flat-files/$(notdir $(3)).tsv --outFileSortedRegions $(3)-sorted-by.bed -bs 10 -p 6 $(4) --beforeRegionStartLength 2000 --afterRegionStartLength 2000"
endef

#}}}

define build-plot-heatmap-command
	mkdir -p $(dir $(2))
	$(eval color := $(call get-colormap,$3))
	${bsub} -n 6 "plotHeatmap -m $1 -out $2 $4 --colorMap '$(color)'"
endef

# Make groups and compute metagene matrices independently # {{{
.PHONY: make-expression-groups
make-expression-groups: ${groups}

${groups}: ${group_by} ${genes_bed}
	awk -F "\t" '{print > "$(strip $(patsubst %.Group-filtered.txt, %, $(@)))."$$2"-filtered.txt"}' $<
#	tr ' ' \\t < ${genes_bed} | bedtools slop -l 2000 -r 0 -s -g ${chromosomes} -i - > $(patsubst %.bed, %-with-promoter.bed, ${genes_bed})
	tr ' ' \\t < ${genes_bed} > $(patsubst %.bed, %-with-tabs.bed, ${genes_bed})
	for f in `ls -d $(patsubst %.Group-filtered.txt, %, $(@))*`; do awk 'FNR==NR { array[$$1]++; next } FS=" " { line = $$10; sub(/^.*=/, "", line); if (line in array) print}' $${f} $(patsubst %.bed, %-with-tabs.bed, ${genes_bed}) > $${f/.txt/.bed}; done

.PHONY: compute-gene-matrix
#compute-gene-matrix: ${genes_bed} ${metagene_matrices} ${metagene_matrices_per_protein}
compute-gene-matrix: make-expression-groups ${metagene_matrices_per_protein} #${metagene_plots}

${deeptools_path}/ses-norm/matrices/%.deeptools-metagene-matrix: $(call KEEP,%,${ses_bigwig}) $(call KEEP,ses-norm,${groups}) .FORCE
	mkdir -p $(@D)
	$(eval x := $(shell ls -d $(patsubst %.Group-filtered.txt, %,$(call KEEP,ses-norm,${groups}))*bed))
	$(call compute-region-matrix,$<,$(strip $(call FILTER_OUT,Group,$(x))),$@,--startLabel 'GS' --endLabel 'GE' --beforeRegionStartLength 2000 --afterRegionStartLength 500)

${deeptools_path}/read-count/matrices/%.deeptools-metagene-matrix: ${deeptools_path}/read-count/%.bw $(call KEEP,read-count,${groups}) 
	mkdir -p $(@D)
	$(eval x := $(shell ls -d $(patsubst %.Group-filtered.txt, %,$(call KEEP,read-count,${groups}))*bed))
	$(call compute-region-matrix,$<,$(strip $(call FILTER_OUT,Group,$(x))),$@,--startLabel 'GS' --endLabel 'GE' --beforeRegionStartLength 2000 --afterRegionStartLength 500)

# }}}

#for p in H3K27ac H3K4me3;
#do tmp=$(realpath paper/groups-main-ones/*.bed | tr "\n" " ");
#~/source/run-bsub.sh -n 8 -M 64000 -R "rusage[mem=30000]" 
#"computeMatrix reference-point -S mm10/bam-compare/read-count/3KO_h0-${p}_pooled_RPGC.bw
#mm10/bam-compare/read-count/3KO_h1-${p}_pooled_RPGC.bw 
#mm10/bam-compare/read-count/3KO_h4-${p}_pooled_RPGC.bw 
#mm10/bam-compare/read-count/3KO_h8-${p}_pooled_RPGC.bw 
#mm10/bam-compare/read-count/3KO_h24-${p}_pooled_RPGC.bw 
#-R ${tmp} --outFileName paper/read-count-main-ones/matrices/${p}_pooled_RPGC.deeptools-NuRD-annotated-matrix
#--outFileSortedRegions paper/read-count-main-ones/matrices/${p}_pooled_RPGC.deeptools-NuRD-annotated-sorted-by.bed -bs 10 -p 8  --beforeRegionStartLength 2000 --afterRegionStartLength 2000"; done
# vim: ft=make
