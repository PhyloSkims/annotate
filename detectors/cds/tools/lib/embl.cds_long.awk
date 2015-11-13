#
# get cds features from embl (long version)
#
# -v FASTA

# @include lib.embl.awk


BEGIN {
  print "#locus locustag genefam gene from to strand nexon length status start stop dnaseq protseq product"

  if (HEADONLY != "") exit(0)

  if (MAXSPAN == "") MAXSPAN = 10000

  if (FASTA == "") Error("No FASTA file specified", 1)
  
  if (! TestPath(FASTA)) Error("Fasta file: '" FASTA "' not found", 1)

  Seq = tolower(ReadFasta(FASTA))
}

/^ID / {
  locus = $2
  gsub(";", "", locus)
  incds = 0
  next
}

/^FT   CDS/ {
  revstrand = match($3, "^complement")
  s = substr($0, 22)
  gsub("^complement", "", s)
  ok = ! match(s, "complement|order")
  nexon = Nexons(s)
  SpanLocation(s, sloc)
  spanlen = sloc[2] - sloc[1] + 1
  len = LenLocation(s)
  ok = ok && (len < MAXSPAN)
  cdsseq = ok ? SeqLocation(Seq, s, revstrand) : "XXX"
  cstart = substr(cdsseq, 1,3)
  cstop  = substr(cdsseq, length(cdsseq)-2)

  gene = "none"
  locustag = "none"
  product = "none"
  translation = "X"
  incds = 1
  next
}

(incds && /^FT   [^ ]/) {
  print locus, locustag, GeneFamily(gene), gene, sloc[1], sloc[2], (revstrand ? "R" : "D"), 
        nexon, len, (ok ? "Ok" : "Error"), cstart, cstop, cdsseq, translation, product
  incds = 0
  next
}

/^FT                   \/gene=/ {
  split($0, a, "=")
  gene = a[2]
  gsub("^[^a-z,A-Z]+", "", gene)
  gsub("\"", "", gene)
  gsub(" ", "_", gene)
  next
}

/^FT                   \/locus_tag=/ {
  split($0, a, "=")
  locustag = a[2]
  gsub("\"", "", locustag)
  gsub(" ", "_", locustag)
  next
}

/^FT                   \/product=/ {
  split($0, a, "=")
  product = a[2]
  gsub("\"", "", product)
  gsub(" ", "_", product)
  next
}

/^FT                   \/translation=/ {
  split($0, a, "=")
  translation = a[2]
  gsub("\"", "", translation)
  gsub(" ", "", translation)
  next
}

END {
  if (incds) {
    print locus, locustag, GeneFamily(gene), gene, sloc[1], sloc[2], (revstrand ? "R" : "D"), 
          nexon, len, (ok ? "Ok" : "Error"), cstart, cstop, cdsseq, translation, product
  }
}
