#!/bin/bash
#
#                   Reverse complement a fasta formated alignment
#
#========================================================================================

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
SCRIPT_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${SCRIPT_DIR}/../../../scripts/bash_init.sh"


function revcomp {
    awk 'function printfasta(seq) {                                  \
            seqlen=length(seq);                                       \
            for (i=1; i <= seqlen; i+=60)                              \
              print substr(seq,i,60);                                 \
         }                                                          \
        function comp(seq) {                                           \
            "echo "seq" | tr acgtACGT tgcaTGCA " | getline res;     \
            close("echo "seq" | tr acgtACGT tgcaTGCA "); \
            return res;                                                \
        }                                                              \
        function rev(seq) {                                            \
            "echo "seq" | rev " | getline res;                         \
            close("echo "seq" | rev ");                                \
            return res;                                                \
        }                                                              \
        function revcomp(seq) {                                        \
            res=rev(comp(seq));                                        \
            return res;                                                \
        }                                                              \
                                                                       \
        (seq) && /^>/ {print head;                                     \
                       printfasta(revcomp(seq));                       \
                       seq=""}                                         \
        /^>/   {head=$0}                                               \
        ! /^>/ {seq=seq$0}                                             \
        END { print head;                                     \
              printfasta(revcomp(seq));                       \
            }' $*
}

revcomp $*