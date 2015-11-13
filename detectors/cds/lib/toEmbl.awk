#
# iff -> embl
#

function FromLoc(i, _local_, s) {
  s = substr(Exon[i]["modif"], 1, 1) "" Exon[i]["from"]
  gsub("=", "", s)
  return s
}

function ToLoc(i, _local_, s) {
  s = substr(Exon[i]["modif"], 2, 1) "" Exon[i]["to"]
  gsub("=", "", s)
  return s
}

function GeneLocation(_local_, d, s) {
  d = Exon[1]["strand"]
  if (d == "+")
      s = FromLoc(1) ".." ToLoc(Nexon)
  else
      s = FromLoc(Nexon) ".." ToLoc(1)
  if (d == "-") s = "complement(" s ")"
  return s
}

function CdsLocation(_local_, d, i, s) {
  d = Exon[1]["strand"]
  if (d == "+") {
    for (i = 1 ; i <= Nexon ; i++)
      s = s "," FromLoc(i) ".." ToLoc(i)
  }
  else {
    for (i = Nexon ; i >= 1 ; i--)
      s = s "," FromLoc(i) ".." ToLoc(i)
  }
  s = substr(s, 2)
  if (Nexon > 1) s = "join(" s ")"
  if (d == "-") s = "complement(" s ")"
  return s
}

function ExonLocation(i, _local_, s) {
  s = FromLoc(i) ".." ToLoc(i)
  if (Exon[i]["strand"] == "-") s = "complement(" s ")"
  return s
}

function Pad(s, len) {
  while (length(s) < len)
    s = s " "
  return s
}

function Feature(feat, loc, _local_, s) {
  s = Pad("FT   " feat, 21)
  print s "" loc
}

function SQualifier(qual, val, _local_, s) {
  s = Pad("FT", 21)
  print s "/" qual "=" val
}

function QQualifier(qual, val) {
  gsub("\"", "''", val)
  SQualifier(qual, "\"" val "\"")
}

function Unk(s) {
  return (s =="" ? "none" : s)
}

#
# rules
#

/^c pass/ {
  PassType = $3
  PassInfo = $NF
  next
}

/^c annot/ {
  GeneName = $4
  Product = $NF
  gsub("_", " ", Product)
  next
}

/^c nclust/ {
  Ngene = $3
  next
}

/^c begin_entry/ {
  Nexon = 0
  delete Exon
  next
}

/^e exon/ {
  Nexon++
  Exon[Nexon]["from"]   = $3
  Exon[Nexon]["to"]     = $4
  Exon[Nexon]["strand"] = $6
  Exon[Nexon]["indels"] = $9 "+" $12
  modif = $15; gsub("\"", "", modif)
  Exon[Nexon]["modif"]  = modif
  next
}

/^e similarity/ {
  split($12, a, "@")
  Simil = a[1] ":" a[2]
  next
}

/^e translate/ {
  Translat = $3
  next
}

/^c end_entry/ {

  GeneName = Unk(GeneName)
  PassType = Unk(PassType)
  
  gname = (Ngene == 1 ? GeneName : GeneName "_" ++Igene)
  locus = ""
  
  Feature("gene", GeneLocation())
  QQualifier("gene", gname)
  QQualifier("locus_tag", locus)
  
  Feature("CDS", CdsLocation())
  SQualifier("codon_start", 1)
  SQualifier("transl_table", 11)
  QQualifier("gene", gname)
  QQualifier("locus_tag", locus)
  QQualifier("product", Product)
  QQualifier("inference", "similar to DNA sequence:" Simil)
  QQualifier("inference", "detect pass:" PassType ":" PassInfo)
  QQualifier("translation", Translat)
  
  if (Nexon > 1) {
    for (i = 1 ; i <= Nexon ; i++) {
      Feature("exon", ExonLocation(i))
      QQualifier("gene", gname)
      QQualifier("locus_tag", locus)
      SQualifier("number", Exon[1]["strand"] == "+" ? i : Nexon-i+1)
    }
  }
}


