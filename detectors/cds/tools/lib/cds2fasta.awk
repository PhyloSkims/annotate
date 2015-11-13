#
# get fasta sequence from cds list 
#
# [-v FIELD=13] CDS sequence
# [-v FIELD=14] Prot Sequence
#

BEGIN {
  if (CHARPERLINE == "") CHARPERLINE = 50
  if (FIELD == "") FIELD = 14
}

/^#/ { next }

{
  name = $1 "@" $2 "@" $3 "@" $5 "@" $6 "@" $7 "@" $8 "@" int($9/3)
  comment = $NF
  PrintFasta($FIELD, name " " comment)
}
