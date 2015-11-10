#
# extend start/stop
#
# -v FASTA START_MODEL STOP_MODEL START_WALK STOP_WALK
#

function UpStart(pos, strand, _local_, seq, i, imax, smax, s) {
  seq = Seq

  if (strand == "-") {
    pos = LenSeq - pos + 1
    seq = RevSeq
  }
  
  imax = 0
  smax = 0
  
  for (i = pos ; i >= Max(1, pos-START_WALK) ; i -= 3) {
     s = substr(seq, i, 3)
     if (s in StopModel) break
     if ((s in StartModel) && (StartModel[s] > smax)) {
       imax = i
       smax = StartModel[s]
     }
  }
  
  if (strand == "-") {
    imax = (imax > 0) ? LenSeq - imax + 1 : imax
  }
  
  return imax
}

#

function DownStop(pos, strand, _local_, seq, i, imax, s) {
  seq = Seq
  
  if (strand == "-") {
    pos = LenSeq - pos + 1
    seq = RevSeq
  }

  imax = 0
  
  for (i = pos ; i < Min(LenSeq, pos+STOP_WALK) ; i += 3) {
     s = substr(seq, i, 3)
     if (s in StopModel) {
       imax = i
       break
     }
  }
  
  if (strand == "-") {
    imax = (imax > 0) ? LenSeq - imax + 1 : imax
  }

  return imax
}

#
# rules
#

BEGIN {

  if (START_MODEL == "") START_MODEL="Models/start.default.frq"
  if (STOP_MODEL == "")  STOP_MODEL="Models/stop.default.frq"
  if (START_WALK == "")  START_WALK=120
  if (STOP_WALK == "")  STOP_WALK=-1
  
  if (! TestPath(FASTA)) Error("Fasta file: '" FASTA "' not found", 1)

  Seq = tolower(ReadFasta(FASTA))
  LenSeq = length(Seq)
  
  RevSeq = RevComplement(Seq)

  if (START_WALK < 0) START_WALK = LenSeq
  if (STOP_WALK < 0) STOP_WALK = LenSeq
  
  if (! TestPath(START_MODEL)) Error("model file: '" START_MODEL "' not found", 2)
  if (! TestPath(STOP_MODEL))  Error("model file: '" STOP_MODEL "' not found", 2)

  ReadModel(START_MODEL, StartModel)
  ReadModel(STOP_MODEL, StopModel)
}

#

/^c begin_entry/ {
  Strand = $7
  nbexon = $9+0
  StartExon = 1     #  always first (even on - strand)
  StopExon = nbexon #  always last  (even on - strand)
  NbExon = 0
}

/^c/ {
  print $0
  next
}

/^e exon/ {
  NbExon++
  if (NbExon == StartExon) {
    pos = (Strand == "+" ? $3+0 : $4+0)
    start = UpStart(pos, Strand)
    mod_start = (start == 0 ? (Strand == "+" ? "<" : ">") : "=")
    start = (start == 0 ? pos : start)
    $(Strand == "+" ? 3 : 4) = start
  } else {
    mod_start = "="
  }
  if (NbExon == StopExon) {
    pos = (Strand == "+" ? $4+0 : $3+0)
    pos += (Strand == "+" ? 1 : -1)
    stop = DownStop(pos, Strand)
    mod_stop = (stop == 0 ? (Strand == "+" ? ">" : "<") : "=")
    last = (stop == 0 ? pos : stop)
    last += (Strand == "+" ? -1 : 1)
    stop = last + (stop == 0 ? 0 : (Strand == "+") ? 3 : -3)
    $(Strand == "+" ? 4 : 3) = stop
  } else {
    mod_stop = "="
  }
  modif = (Strand == "+" ? mod_start "" mod_stop : mod_stop "" mod_start)
  print $0, "; modifier \"" modif "\""
  next
}

/^e (intron|splice|similarity)/ {
  print $0
  next
}
