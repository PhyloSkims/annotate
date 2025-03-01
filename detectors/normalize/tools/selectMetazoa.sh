#!/bin/bash

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${THIS_DIR}/../../../scripts/bash_init.sh"


(                                          \
	for f in $1/*.gbk* ; do                       \
		if [[ "$f" =~ \.gz$ ]] ; then      \
			GREP=zgrep;                    \
	    else                               \
			GREP=grep;                     \
	    fi;                                \
		${GREP} -H -A 1 '  ORGANISM' $f;   \
    done                                   \
) | \
  grep -B 1 Metazoa | \
  $AwkCmd '{print $1}' | \
  grep '\.gbk' | \
  sed -E 's/(^.*\.gbk(.gz)?).$/\1/' | \
  uniq