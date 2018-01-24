#
# translate CDSs from iff file
#

BEGIN {
  Seq = ReadFasta(FASTA)
}

/^c end_entry/ {
  if (RevStrand) Cds = RevComplement(Cds)
  DNA = substr(Cds, 1, length(Cds)-3)
  print ">" genename 
  for (i=1; i<=length(DNA); i+=60)
  	print substr(DNA,i,60)
}

/^c annot/ {
	genename=$3
}

/^c begin_entry/ {
  Cds = ""
  Iexon = 0
  next
}

/^e exon/ {
  RevStrand = ($6 == "-")
  if (++Iexon == 1) { # first is exon with start (even on - strand)
    Modif = $15
    gsub("\"", "", Modif)
    Modif = (RevStrand ? substr(Modif, 2, 1) : substr(Modif, 1, 1))
  }
  if (RevStrand)
    Cds = SubSeq(Seq, $3, $4) "" Cds
  else
    Cds = Cds "" SubSeq(Seq, $3, $4)
  next
}
