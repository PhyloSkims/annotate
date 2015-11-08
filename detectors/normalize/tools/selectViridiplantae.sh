#!/bin/bash

grep -A 1 '  ORGANISM' $* | \
  grep -B 1 Viridiplantae | \
  gawk '{print $1}' | \
  grep '\.gbk' | \
  sed -E 's/(^.*\.gbk).$/\1/' | \
  uniq