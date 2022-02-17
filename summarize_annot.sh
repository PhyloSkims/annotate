#!/bin/bash

egrep "^FT   (CDS|rRNA|tRNA) |/gene=" $1 \
  | awk '
        /\/gene=/ && (gene!=1) {
            print $0; 
            gene=1
            } 
        /^FT   (CDS|rRNA|tRNA)/ {
            print $0;
            gene=0
            }' 