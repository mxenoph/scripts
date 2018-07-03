#!/usr/bin/env bash

# https://gist.github.com/klmr/e51822f624292681ba98
# Inspired by <http://blog.melski.net/2010/11/30/makefile-hacks-print-the-value-of-any-variable/>

usage() {
    echo >&2 "Usage: $0 [-f makefile] variables..."
    echo >&2
    echo >&2 "Optional arguments:"
    echo >&2 "  -f makefile:  path to the makefile"
    echo >&2
    echo >&2 "Positional arguments:"
    echo >&2 "  variables...: one or more names of variables to print"
}

printenv_make=$(cat <<'MAKE'
debugprint-%:
	@echo '$*=$($*)'
	@echo '  origin = $(origin $*)'
	@echo '  flavor = $(flavor $*)'
	@echo '   value = $(value  $*)'
MAKE
)

if [ "$1" == "-f" ]; then
    shift
    filename="$1"
    shift
else
    filename=""
    if [ -f GNUMakeFile ]; then
        filename=GNUMakeFile
    elif [ -f makefile ]; then
        filename=makefile
    elif [ -f Makefile ]; then
        filename=Makefile;
    fi

    if [ ! -n "$filename" ]; then
        echo >&2 "No makefile found"
        exit 1
    fi
fi

if [ $# -eq 0 ]; then
    usage
    exit 1
fi

vars=()
for var in "$@"; do
    vars+=("debugprint-$var")
done

make -f <(echo "$printenv_make") -f "$filename" "${vars[@]}"
