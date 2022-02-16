#!/bin/bash
#

# -- CAUTION -- Works as long than the script 
#               is not called through a symlink
THIS_DIR="$(dirname ${BASH_SOURCE[0]})"
source "${THIS_DIR}/../../../scripts/bash_init.sh"


  $AwkCmd 'function printfasta(seq) {                                            \
             seqlen=length(seq);                                             \
             for (i=1; i <= seqlen; i+=60)                                   \
                 print substr(seq,i,60);                                     \
             }                                                               \
       function comp(seq) {                   						                   \
              "echo "seq" | tr acgtACGT tgcaTGCA " | getline res; 	         \
              close("echo "seq" | tr acgtACGT tgcaTGCA ");                   \
              return res;                        						                 \
             }                                      						             \
       function rev(seq) {                    						                   \
              "echo "seq" | rev " | getline res; 						                 \
              close("echo "seq" | rev ");                                    \
              return res;                        						                 \
             }                                      						             \
       function revcomp(seq) {                						                   \
              res=rev(comp(seq));                						                 \
              return res;                        						                 \
             }                                      						             \
                                                                             \
       /^LOCUS / {AC=$2; sequence=""; seqon=0; FROM="";TO=""}                \
       /^     misc_feature/ {LOCUS=$2; STRAND=1}                             \
       /^     misc_feature/ && /complement/ {STRAND=0;                       \
                                             sub("complement\\(","",LOCUS);  \
                                             sub("\\)","",LOCUS);    \
                                            }                      \
       /large single copy/ && /LSC/ {split(LOCUS,POS,".");         \
                                     FROM=POS[1];                  \
                                     TO=POS[3];                    \
                                     LENGTH=TO-FROM+1              \
                                    }                              \
       /^ORIGIN/ {seqon=1}                                         \
       /^ *[1-9][0-9]* [a-z ]+$/ && seqon {seq=$2 $3 $4 $5 $6 $7;  \
                                           gsub("[^acgt]","n",seq);\
                                           sequence=sequence seq   \
                                          }                        \
       /^\/\// && FROM && (LENGTH > 60000) && (LENGTH < 100000)    \
                        {print ">LSC_"AC" Strand="STRAND";",       \
                               "cut="FROM".."TO";",                \
                               "seq_length="LENGTH";";             \
                         SS=substr(sequence,FROM,LENGTH);          \
                         if (! STRAND)                             \
                           SS=revcomp(SS);                         \
                         printfasta(SS);                           \
       }             \
      ' $*
