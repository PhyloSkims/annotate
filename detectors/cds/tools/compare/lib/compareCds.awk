#
#
#

function Min(a, b) {
  return (a < b ? a : b)
}

function Max(a, b) {
  return (a > b ? a : b)
}

function Align(s1, s2, _local_, d, l) {
  if (s1 == s2) return 100
  d = AlignNWS(s1, s2, Identity)
  l = Max(length(s1), length(s2))
  return int((l - d) * 100 / l)
}

BEGIN {
  PROCINFO["sorted_in"] = "@ind_num_asc"
  IdentityMatrix("ABCDEFGHIJKLMNOPQRSTUVWXYZ*", Identity)
}

BEGINFILE {
  NFile++
  File[NFile] = FILENAME
}

/^#/ { next }

{
  strand = $5
  stop = (strand == "D" ? $4 : $3)
  Stop[stop]++
  i = ++NRec[NFile]
  Rec[NFile][i]["record"] = $0
  Rec[NFile][i]["genefam"] = $1
  Rec[NFile][i]["gene"] = $2
  Rec[NFile][i]["from"] = $3
  Rec[NFile][i]["to"] = $4
  Rec[NFile][i]["strand"] = $5
  Rec[NFile][i]["nexon"] = $6
  Rec[NFile][i]["length"] = $7
  Rec[NFile][i]["protseq"] = $9
  if (NFile == 1)
    Indx1[stop] = i
  else
    Indx2[stop] = i
}

END {
  for (st in Stop) {
    if (Indx1[st])
      print "FILE1 " Rec[1][Indx1[st]]["record"]
    else 
      print "FILE1 NONE"

    if (Indx2[st])
      print "FILE2 " Rec[2][Indx2[st]]["record"]
    else 
      print "FILE2 NONE"
      
    if (Indx1[st] && Indx2[st]) {
      fm = Rec[1][Indx1[st]]["genefam"]
      id = Align(Rec[1][Indx1[st]]["protseq"], Rec[2][Indx2[st]]["protseq"])
      printf("MATCH %s ID %d ", fm, id)
      if (id == 100)
        status = "CORRECT"
      else if (id >= 90)
        status = "ALMOST_CORRECT"
      else if (id >= 80)
        status = "ACCEPTABLE"
      else
        status = "WRONG"
      if (status != "CORRECT") {
        if (Rec[1][Indx1[st]]["nexon"] != Rec[2][Indx2[st]]["nexon"])
          status = status ".BAD_NBEXON"
        start1 = Rec[1][Indx1[st]]["strand"] == "D" ? Rec[1][Indx1[st]]["from"] : Rec[1][Indx1[st]]["to"]
        start2 = Rec[2][Indx2[st]]["strand"] == "D" ? Rec[2][Indx2[st]]["from"] : Rec[2][Indx2[st]]["to"]
        if (start1 != start2)
          status = status ".BAD_START"
        else
          status = status ".BAD_JUNCTION"
      }
      print status
    }
    else if (Indx1[st]) {
      fm = Rec[1][Indx1[st]]["genefam"]
      print "MATCH " fm " ID 0 MISSED.WRONG_STOP"
    }
    else if (Indx2[st]) {
      fm = Rec[2][Indx2[st]]["genefam"]
      print "MATCH " fm " ID 0 OVERPRED.WRONG_STOP"
    }
    print ""
  }
}
