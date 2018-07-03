.DEFAULT_GOAL := show-help
.DELETE_ON_ERROR:

# See <https://gist.github.com/mschubert/a0e4f3aeaf3558431890> for explanation.
.PHONY: .FORCE
.FORCE: ;

define EXT_DEP
$1/$2: .FORCE
	make -f $(1).make $(2)
endef

ifeq (n,$(findstring n,$(firstword -$(MAKEFLAGS))))
	ext_dep = $(NOOP)
else
	ext_dep = $(foreach i,$(2),$(eval $(call EXT_DEP,$(1),$(i))))
endif

# vim: ft=make
