#
# select best cluster(s) 
#
# -v MAX_SPAN ALLOW_STOP EXCLUDE
#
# 
#

BEGIN {
  PROCINFO["sorted_in"] = "@ind_num_asc"
  if (MAX_SPAN == "") MAX_SPAN = 10000
  if (ALLOW_STOP == "") ALLOW_STOP  = 0
  if (EXCLUDE == "") EXCLUDE = 0
}

/^# --- START OF GFF DUMP ---/ {
  State = "gff"
  NbEntry++
  next
}

/^# --- END OF GFF DUMP ---/ {
  State = 0
  next
}

/^C4 Alignment:/ {
  for (i = 1 ; i <= 8 ; i++)
    getline
  State = "align"
  Align = ""
  next
}

/^#/ { next }

(State == 0) { next }

(State == "gff") && ($3 == "gene") {
  span = $5 - $4 + 1
  valid =    (span < MAX_SPAN) \
          && (ALLOW_STOP || (Align !~ "\\*")) \
          && ((EXCLUDE == 0) || (! (Organism ~ EXCLUDE)))
  Entry[NbEntry]["valid"] =  valid
  Entry[NbEntry]["from"] = $4+0
  Entry[NbEntry]["to"] = $5+0
  Entry[NbEntry]["score"] = $6+0
  Entry[NbEntry]["strand"] = $7
  Entry[NbEntry]["stop"] = (Align !~ "\\*")
  if (valid) {
    for (i = $4+0 ; i <= $5+0; i++)
      Cover[i] = 1
  }
}

(State == "gff") && ($3 == "exon") {
  Entry[NbEntry]["nbexon"]++
}

(State == "gff") {
  n = ++Entry[NbEntry]["nbline"]
  $1=$2=""
  Entry[NbEntry][n] = $0
}

(State == "align") && /^vulgar/ {
  State = 0
  split($2, a, "@")
  Organism = a[1]
  next
}

(State == "align") {
  getline; getline
  Align = Align "" $0
  getline; getline
  next
}

END {
  ClustScoreMax=0
  # make clusters
  pi = -1
  for (i in Cover) {
    if (i+0 > pi+1)
      Clust[++NbClust]["from"] = i
    pi = Clust[NbClust]["to"] = i
  }
  
  # get highest scoring clusters
  for (i = 1 ; i <= NbEntry ; i++) {
    valid = Entry[i]["valid"]
    if (! valid) continue
    clusno = 0
    for (j = 1; j <= NbClust; j++) {
      if ((Entry[i]["from"] >= Clust[j]["from"]) && (Entry[i]["to"] <= Clust[j]["to"]))
        clusno = j
    }
    valid = (clusno != 0)
    if (! valid) continue
    
    score = Entry[i]["score"]
    if (score > Clust[clusno]["score"]+0) {
      Clust[clusno]["score"] = score
      Clust[clusno]["strand"] = Entry[i]["strand"]
      Clust[clusno]["entry"] = i
      
      if (score > ClustScoreMax)
      	ClustScoreMax=score
    }
  }

  
  # Consider only cluster with at least 95% of the best score
  
  NbClustOk=0
  for (i = 1 ; i <= NbClust ; i++) {
    if (Clust[i]["score"] >= ClustScoreMax * 0.95) 
    	NbClustOk++
    else
    	Clust[i]["score"]=0

  }


  # print cluster info
  print "c nclust", NbClustOk+0
  for (i = 1 ; i <= NbClust ; i++) {
  	if (Clust[i]["score"] == 0) continue
    print "c cluster", i, "from", Clust[i]["from"], "to", Clust[i]["to"],\
          "strand", Clust[i]["strand"], "score", Clust[i]["score"]
  }
  
  # print best clusters
  for (i = 1 ; i <= NbClust ; i++) {
    if (Clust[i]["score"] == 0) continue
    j = Clust[i]["entry"]
    s = Clust[i]["strand"]
    n = Entry[j]["nbline"]
    ne = Entry[j]["nbexon"]
    print "c begin_entry", j, "cluster", i, "strand", s, "nbexon", ne
    for (k = 1 ; k <= n ; k++) {
      entry = Entry[j][k]
      gsub("^ +", "", entry)
      print "e", entry
    }
    print "c end_entry", j, "cluster", i
  }
}


