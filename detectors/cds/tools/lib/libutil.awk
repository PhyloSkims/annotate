#
# utilities library
#

END {
  if (_ForceExit_) exit(_ForceExit_)
}

function Notify(key, msg) {
    print "# " key " " msg >> "/dev/stderr"
}

function Info(msg) {
    Notify("info", msg)
    return 1
}

function Warning(msg) {
    Notify("warning", msg)
    return 0
}

function Exit(status) {
    exit(_ForceExit_ = status)
}

function Error(msg, status) {
    Notify("error", msg)
    Exit(status ? status : 1)
}

function Assert(condition, msg, status) {
    if (! condition) {
      msg = FILENAME ":" FNR ": " msg
      return status ? Error(msg, status) : Warning(msg)
    }
}

function Min(a, b) {
  return (a < b ? a : b)
}

function Max(a, b) {
  return (a > b ? a : b)
}

function Trim(s, regexp) {
  if (regexp == 0) regexp = "[ \t]+"
  gsub("^" regexp "|" regexp "$", "", s)
  return s
}

function TestPath(path, test, _local_, stat) {
  if (test == 0) test = "-f"
  if (Trim(path) == "")
    return 0 # because of a bug in 'test'

  gsub(" ", "\\ ", path) # protect spaces
  stat = system("test " test " " path)
  return stat ? 0 : 1
}

# ---

function Strip(s) {
  gsub("complement|join|order|\\)|\\(|<|>| ", "", s)
  return s
}

function GetLoc(s, loc, _local_, a, tmp) {
  delete loc
  loc[1] = loc[2] = 0
  split(s, loc, "\\.\\.")
  if (loc[1] > loc[2]) {
    tmp = loc[1]
    loc[1] = loc[2]
    loc[2] = tmp
  }
}

function ParseLocation(s, locs, _local_, i, na, a, loc) {
  delete locs
  s = Strip(s)
  na = split(s, a, ",")
  for (i = 1 ; i <= na ; i++) {
    GetLoc(a[i], loc)
    locs[i][1] = loc[1]
    locs[i][2] = loc[2]
  }
  return na
}

function SpanLocation(s, sloc, _local_, i, na, locs) {
  delete sloc
  na = ParseLocation(s, locs)
  sloc[1] = (na > 0 ? locs[1][1] : 0)
  sloc[2] = (na > 0 ? locs[1][2] : 0)
  for (i = 2 ; i <= na ; i++) {
    sloc[1] = Min(sloc[1], locs[i][1])
    sloc[2] = Max(sloc[2], locs[i][2])
  }
}

function LenLocation(s, _local_, i, na, locs, len) {
  len = 0
  na = ParseLocation(s, locs)
  for (i = 1 ; i <= na ; i++) {
    len += locs[i][2] - locs[i][1] + 1
  }
  return len
}

function Nexons(s, _local_, a) {
  s = Strip(s)
  return split(s, a, ",")
}

function Reverse(s, _local_, i, n, rs) {
  rs = "";
  n = length(s);
  for (i = n ; i >= 1 ; i--)
    rs = rs "" substr(s, i, 1)
  return rs;
}

function RevComplement(seq, _local_, n, i, c, rs) {
  n = length(seq)
  rs = ""
  for (i = n ; i >= 1 ; i--) {
    c = substr(seq, i, 1)
    rs = rs "" (_DnaC[c] ? _DnaC[c] : "X")
  }
  return rs
}

function SubSeq(seq, from, to, revstrand) {
  seq = substr(seq, from, to-from+1)
  if (revstrand) seq = RevComplement(seq)
  return seq
}

function SeqLocation(seq, s, revstrand, _local_, sloc, i, na, locs) {
  sloc = ""
  na = ParseLocation(s, locs)
  for (i = 1 ; i <= na ; i++) {
    sloc = sloc "" SubSeq(seq, locs[i][1], locs[i][2], 0)
  }
  return (revstrand ? RevComplement(sloc) : sloc)
}

function GeneFamily(s) {
  s = tolower(s)
  gsub("(_|-)[0-9]+$", "", s)
  gsub("(_|-)(a|b|c|i)+$", "", s)
  gsub("'+$", "", s)
  gsub("/.+$", "", s)
  if (match(s, "[^a-z,0-9]")) s = "none"
  return s
}

# ---

function ReadFasta(file, _local_, seq, context) {
  context = $0
  seq = ""
  getline < file
  while(getline < file) seq = seq "" $0
  $0 = context
  return seq    
}

function PrintFasta(seq, name, CharPerLine, _local_, i, len, rest, chunk) {
  if (CharPerLine <= 0) CharPerLine = 60
  print ">" name
  len = length(seq)
  rest = len % CharPerLine
  chunk  = (len - rest) / CharPerLine
  for (i = 0 ; i < chunk ; i++) 
    print substr(seq, (i*CharPerLine) + 1, CharPerLine)
  if (rest != 0)
    print substr(seq, (chunk*CharPerLine) + 1, rest)
}


#
# constants
#

BEGIN {
  # complementary of _IupacDna
  _DnaC["A"] = "T"; _DnaC["C"] = "G"; _DnaC["G"] = "C"; _DnaC["T"] = "A"
  _DnaC["R"] = "Y"; _DnaC["Y"] = "R"; _DnaC["M"] = "K"; _DnaC["K"] = "M"
  _DnaC["W"] = "W"; _DnaC["S"] = "S"; _DnaC["B"] = "V"; _DnaC["V"] = "B"
  _DnaC["D"] = "H"; _DnaC["H"] = "D"; _DnaC["N"] = "N"; _DnaC["X"] = "X"
  _DnaC["a"] = "t"; _DnaC["c"] = "g"; _DnaC["g"] = "c"; _DnaC["t"] = "a"
  _DnaC["r"] = "y"; _DnaC["y"] = "r"; _DnaC["m"] = "k"; _DnaC["k"] = "m"
  _DnaC["w"] = "w"; _DnaC["s"] = "s"; _DnaC["b"] = "v"; _DnaC["v"] = "b"
  _DnaC["d"] = "h"; _DnaC["h"] = "d"; _DnaC["n"] = "n"; _DnaC["x"] = "x"
}
