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
	FrameShift=0
	delete Exon
	next
}

/^e exon/ {
  Nexon++
  Exon[Nexon]["from"]   = $3
  Exon[Nexon]["to"]     = $4
  Exon[Nexon]["strand"] = $6

	match($0, / insertions +([0-9]+) +/, arr)
	if (arr[1])
		insertions=arr[1]
	else
		insertions=0
		
	match($0, / insertions +([0-9]+) +/, arr)
	if (arr[1])
		deletions=arr[1]
	else
		deletions=0
				
	Exon[Nexon]["indels"] = insertions "+" $deletions
		
	match($0, / modifier +"([^"]*)"/, arr)
	if (arr[1])
		modif=substr(arr[1],1,1)
	else
		modif=""
	if (modif == "=")
		modif=""
		
	Exon[Nexon]["modif"]  = modif

	match($0, / frameshifts +([-0-9]+)/, arr)
	if (arr[1]) {
		FrameShift=FrameShift+1
		Exon[Nexon]["frameshift"] = arr[1]
	}
	else
		Exon[Nexon]["frameshift"] = 0
	
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
	  
	# Feature("gene", GeneLocation())
	#QQualifier("gene", gname)
	#QQualifier("locus_tag", locus)
	#if (FrameShift)
	#	QQualifier("pseudogene","unknown")
		  
	Feature("CDS", CdsLocation())
	SQualifier("codon_start", 1)
	SQualifier("transl_table", 11)
	QQualifier("gene", gname)
	QQualifier("locus_tag", locus)
	if (FrameShift) {
		QQualifier("pseudogene","unknown")
		if (FrameShift > 1)
			QQualifier("note","nonfunctional due to frameshifts in "FrameShift" exons")
		else
			QQualifier("note","nonfunctional due to a frameshift")
		}
	QQualifier("product", Product)
	QQualifier("inference", "similar to DNA sequence:" Simil)
	# QQualifier("inference", "org.annot -- detect pass:" PassType ":" PassInfo)
	if (FrameShift==0) {
		if (match(Translat,/\*/)>0) {
			QQualifier("pseudogene","unknown")
			QQualifier("note","nonfunctional due to stop codon")
		}
		
		QQualifier("translation", Translat)
	}
		
	  
	if (Nexon > 1) {
		for (i = 1 ; i <= Nexon ; i++) {
	    	Feature("exon", ExonLocation(i))
	    	QQualifier("gene", gname)
	    	QQualifier("locus_tag", locus)
	    	
	    	if (Exon[i]["frameshift"]) {
				QQualifier("pseudogene","unknown")
				if (Exon[i]["frameshift"] > 0)
					QQualifier("note","frameshifted by insertion of " Exon[i]["frameshift"] " bp")
				else
					QQualifier("note","frameshifted by deletion of " -Exon[i]["frameshift"] " bp")
			}
	    			    		
	    	SQualifier("number", Exon[1]["strand"] == "+" ? i : Nexon-i+1)
	    }
	  }
}


