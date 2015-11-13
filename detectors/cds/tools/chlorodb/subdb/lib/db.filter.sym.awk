#
#
#

function Check(seq) {
  if (seq == "") return 0
  gsub("[ACDEFGHIKLMNPQRSTVWXY\n]+", "", seq)
  return (length(seq) == 0)
}

/^>/ {
  if (Check(Seq)) {
    print Name
    printf("%s", Seq)
  }
  Name = $0
  Seq = ""
  next
}

{
  Seq = Seq "" $0 "\n"
}

END {
  if (Check(Seq)) {
    print Name
    printf("%s", Seq)
  }
}
