#
# translate CDSs from iff file
#

BEGIN {
  Seq = ReadFasta(FASTA)
}

/^c end_entry/ {
  if (RevStrand) Cds = RevComplement(Cds)
  Prot = Translate(substr(Cds, 1, length(Cds)-3),Modif,"chloroplast")

  if (Modif == "=" && substr(Prot, 1,1) != "M") {
    Prot = "M" substr(Prot, 2,length(Prot))
    if (CdsStartPosFull==1) {
      print "e trans_exception "CdsStartPos" Met"
    } else {
      print "e trans_exception " "ERROR" " Met"
    }
  }
  print "e translate " Prot
}

{
  print $0
}

/^c begin_entry/ {
  Cds = ""
  Iexon = 0
  CdsStartPosFullv = 0
  next
}

/^e exon/ {
  print "d cds " CDS
  RevStrand = ($6 == "-")
  if (++Iexon == 1) { # first is exon with start (even on - strand)
    Modif = $15
    gsub("\"", "", Modif)
    Modif = (RevStrand ? substr(Modif, 2, 1) : substr(Modif, 1, 1))
  }
  if (RevStrand){
    Cds = SubSeq(Seq, $3, $4) "" Cds
    if (Iexon==1) {
      if(($4 - $3 +1) >= 3) {
        CdsStartPos="complement("$4-2".."$4")"
        CdsStartPosFull=1
      }
    }}
  else{
    Cds = Cds "" SubSeq(Seq, $3, $4)
    if (Iexon==1) {
      if (Iexon==1) {
        if(($4 - $3 +1) >= 3) {
          CdsStartPos=$3".."$3+2
          CdsStartPosFull=1
        }
    }
    }  }
  next
}
