#
#

function visit(u, i, _local_, v) {
  if (Visited[u]) return
  Visited[u] = i
  for (v in Edge[u]) {
    visit(v, i)
  }
}

/^#/ { next }

{
  Node[$1]++
  Node[$2]++
  Edge[$1][$2]++
  Edge[$2][$1]++
}

END {

  for (u in Node) {
    if (Visited[u]) continue
    visit(u, ++NbComp)
  }

  for (u in Node)
    print u, Visited[u]
  
}