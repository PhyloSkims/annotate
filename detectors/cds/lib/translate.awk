#
# translate CDSs from iff file
#

BEGIN {
  Seq = ReadFasta(FASTA)
}

/^c end_entry/ {
  if (RevStrand) Cds = RevComplement(Cds)
  Prot = Translate(substr(Cds, 1, length(Cds)-3))
  if (Modif == "=")
    Prot = "M" substr(Prot, 2)
  print "e translate " Prot
}

{
  print $0
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
