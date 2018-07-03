# Applications used
bsub = /homes/mxenoph/source/run-bsub.sh -K -q research-rh7
java = /ebi/research/software/Linux_x86_64/opt/java/jdk1.6/bin/java -Xmx2g -jar
picard = /nfs/research2/bertone/software/picard-tools-1.76

# Other parameters
memlimit = 64000

#.PHONY: FORCE
#FORCE: ;

# Fast way to force a rule for debugging
#.PHONY: .FORCE
#.FORCE: ;

# Rule to print any variable in this makefile by invoking print-VARNAME
print-%:
	@echo $* = $($*)

# Print make_version
version:
	@echo GNU Make v$(MAKE_VERSION)

# macros: foreach element in $(2) i.e. mapped_reads# {{{
# search for string $(1) i.e. nput and if true replace element
# with empty
FILTER_OUT = $(foreach element, $(2),\
			 $(if $(findstring $(1),$(element)),,$(element)))
KEEP = $(foreach element, $(2),\
	   $(if $(findstring $(1),$(element)),$(element)))
# }}}

#Showing help# {{{

.PHONY: show-help
show-help:
	@echo "$$(tput bold)Available rules:$$(tput sgr0)"
	@echo
	@sed -n "/^##/ { \
		h; \
		n; \
		s/:.*//; \
		G; \
		s/^/$$(tput setaf 6)/; \
		s/\\n##/$$(tput sgr0)---/; \
		p; \
	}" ${MAKEFILE_LIST} \
	| sort --ignore-case \
	| awk 'BEGIN {FS = "---"} { printf "%-30s\t%s\n", $$1, $$2 }'
# }}}

# Macros specific to peak calling for getting appropriate control # {{{

get-base = $(shell echo $(notdir $1) | sed 's/-.*/-/g')
get-ctrl-type = $(shell echo $(call KEEP,$(shell echo $(notdir $1) | sed 's/.*-//g' | sed 's/_.*//'),\
	${dict_chip_ctrl_type}) | sed 's/.*-//g')
get-ctrl = $(call KEEP,$(3),$(filter $(1)/$(2)%, $(4)))
format-ctrl = $(shell if [ -z $(strip $(1)) ] ; then echo ''; else echo '-c $(1)'; fi)
#}}}

# vim:ft=make
