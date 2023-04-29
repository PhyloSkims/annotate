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

#

function Min(a, b) {
  return (a < b ? a : b)
}

function Max(a, b) {
  return (a > b ? a : b)
}

#

function Trim(s, regexp) {
  if (regexp == 0) regexp = "[ \t]+"
  gsub("^" regexp "|" regexp "$", "", s)
  return s
}

function ShieldPath(path) {
  gsub(" ", "\\ ", path)
  return path
}

function TestPath(path, test, _local_, stat) {
  if (test == 0) test = "-f"
  if (Trim(path) == "")
    return 0 # because of a bug in 'test'
  stat = system("test " test " " ShieldPath(path))
  return stat ? 0 : 1
}

#

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

#

function _AssertCode(name, _local_, n, i1, i2, i3, b1, b2, b3) {
  if (_InitCod[name] != 0) return 1
  for (i1 = 1 ; i1 <= 4 ; i1++) {
    b1 = substr(_NucOrder, i1, 1)
    for (i2 = 1 ; i2 <= 4 ; i2++) {
      b2 = substr(_NucOrder, i2, 1)
      for (i3 = 1 ; i3 <= 4 ; i3++) {
        b3 = substr(_NucOrder, i3, 1)
        _GenCod[name][b1 ""  b2 "" b3] = substr(_Cod2Aa[name], ++n, 1)
        if (name in  _Cod2Start) {
          _GenStart[name][b1 ""  b2 "" b3] = substr(_Cod2Start[name], n, 1)
          if (_GenStart[name][b1 ""  b2 "" b3] == "-") {
            _GenStart[name][b1 ""  b2 "" b3] = _GenCod[name][b1 ""  b2 "" b3]
          }
        } else {
          _GenStart[name][b1 ""  b2 "" b3] = _GenCod[name][b1 ""  b2 "" b3]
        }
      }
    }
  }
  _InitCod[name] = 1
  return 1
}

function Translate(seq, modif, codname, _local_, n, i, c, p) {
  if (codname == 0) codname = "universal"
  print "u codname " codname " " _Cod2Aa[name]
  p= ""
  _AssertCode(codname)
  seq = toupper(seq)

  n = length(seq)
  for (i = 1 ; i <= n ; i += 3) {
    c = substr(seq, i, 3)
    if (i == 1 && modif == "=" )
       aa = ((c in _GenStart[codname]) ? _GenStart[codname][c] : "X")
    else 
       aa = ((c in _GenCod[codname]) ? _GenCod[codname][c] : "X")

    p = p "" aa
  }

  return p
}

#

function SubSeq(seq, from, to, revstrand) {
  seq = substr(seq, from, to-from+1)
  if (revstrand) seq = RevComplement(seq)
  return seq
}

#

function ReadFasta(file, _local_, seq, context) {
  context = $0
  seq = ""
  getline < file
  while(getline < file) seq = seq "" $0
  $0 = context
  return seq    
}

#

function ReadModel(file, a, _local_, line, context) {
  context = $0
  delete a
  while(getline < file)
    if (! ($0 ~ "^#")) a[$1] = $2
  $0 = context
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

  # genetic codes
  _NucOrder = "ACGT"
                                #1  AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT
                                #2  AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
                                #3  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
  _Cod2Aa["universal"]           = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"
  _Cod2Aa["chloroplast"]         = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"
  _Cod2Aa["mito-yeast"]          = "KNKNTTTTRSRSMIMIQHQHPPPPRRRRTTTTEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"
  _Cod2Aa["mito-vertebrates"]    = "KNKNTTTT*S*SMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"
  _Cod2Aa["mito-insects"]        = "KNKNTTTTSSSSMIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"
  _Cod2Aa["mito-echinoderms"]    = "NNKNTTTTSSSSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"  
  _Cod2Aa["mycoplasma"]          = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSWCWCLFLF"
  _Cod2Aa["ciliata"]             = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVVQYQYSSSS*CWCLFLF"
  _Cod2Aa["euplotes"]            = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRLLLLEDEDAAAAGGGGVVVV*Y*YSSSSCCWCLFLF"
  _Cod2Aa["candida"]             = "KNKNTTTTRSRSIIMIQHQHPPPPRRRRXLSLEDEDAAAAGGGGVVVV*Y*YSSSS*CWCLFLF"

                                #1  AAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTT
                                #2  AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTT
                                #3  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
  _Cod2Start["universal"]        = "------------M-----------------M-------------------------------M-"
  _Cod2Start["chloroplast"]      = "------------MMMM--------------M---------------M---------------M-"
}
