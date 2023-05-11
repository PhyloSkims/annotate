#!/bin/bash

resume_ft() {
    local input=$1
    grep -E "^FT   (CDS|rRNA|tRNA) |/gene=" ${input} \
    | awk '
            /\/gene=/ && (gene!=1) {
                print $0; 
                gene=1
                } 
            /^FT   (CDS|rRNA|tRNA)/ {
                print $0;
                gene=0
                }' 
}

count_tRNA() {
    local input=$1

    resume_ft $input \
    | grep '^FT   tRNA' \
    | wc -l
}




printf "tRNA genes : %d\n" "$(count_tRNA $1)"

