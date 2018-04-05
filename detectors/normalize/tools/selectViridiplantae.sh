#!/bin/bash


(                                          \
	for f in $* ; do                       \
		if [[ "$f" =~ \.gz$ ]] ; then      \
			GREP=zgrep;                    \
	    else                               \
			GREP=grep;                     \
	    fi;                                \
		${GREP} -H -A 1 '  ORGANISM' $f;   \
    done                                   \
) | \
  grep -B 1 Viridiplantae | \
  gawk '{print $1}' | \
  grep '\.gbk' | \
  sed -E 's/(^.*\.gbk(.gz)?).$/\1/' | \
  uniq