#!/bin/bash
#
#                           splitgbk.sh:
#                           Split a gbk file in multiple files
#                           each containing a single sequence
#
#========================================================================================

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${THIS_DIR}/../../../scripts/bash_init.sh"

inputfile=$1
dest=${inputfile/.*/}

mkdir -p $dest

$AwkCmd -v dest="$dest" '/^LOCUS/ {
            AC=$2;
            destfile = sprintf("%s/%s.gbk", dest, AC);
        }
        { 
            print $0 >> destfile
        }
        /^\/\// {
            close(destfile);
        }
        ' $inputfile