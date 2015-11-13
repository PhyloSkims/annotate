#
#

function getOrg(s, _local_, a, na, org) {
  na = split(s, a, "/")
  na = split(a[na], a, "\\.")
  return a[1]  
}


BEGIN {
  PROCINFO["sorted_in"] = "@ind_num_asc"
  print "#org tot cor alcor acc wrong over misstot misschlo missoth"
}

/MISSED in ChloroDB/ {
  org = getOrg($1)
  Org[org]++
  Cnt[org]["MISSCHLORO"] = $2
  next
}

/MISSED not in ChloroDB/ {
  org = getOrg($1)
  Org[org]++
  Cnt[org]["MISSNOTCHLORO"] = $2
  next
}

/^#/ { next }

/^.*:MATCH/ {
  org = getOrg($1)
  Org[org]++
  split($NF, a, "\\.")
  Cnt[org][a[1]]++
}

END {
  for (org in Org) {
    Cnt[org]["TOTAL"] =   Cnt[org]["CORRECT"] + Cnt[org]["ALMOST_CORRECT"] \
                        + Cnt[org]["ACCEPTABLE"] + Cnt[org]["WRONG"]       \
                        + Cnt[org]["MISSED"]
  }
  for (org in Org) {
    print org, Cnt[org]["TOTAL"]+0, Cnt[org]["CORRECT"]+0,          \
          Cnt[org]["ALMOST_CORRECT"]+0, Cnt[org]["ACCEPTABLE"]+0,   \
          Cnt[org]["WRONG"]+0, Cnt[org]["OVERPRED"]+0,              \
          Cnt[org]["MISSED"]+0,                                     \
          Cnt[org]["MISSCHLORO"]+0, Cnt[org]["MISSNOTCHLORO"]+0
          
  }

}
