#

/^>/ {
  split($1, a, "@")
  ok = a[3] ~ PAT 
}

ok {
  print $0
}
