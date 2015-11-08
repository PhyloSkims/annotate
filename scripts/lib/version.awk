#
# check version number
#
# -v NUM=<curversion>
# -v REF=<minversion>
#

BEGIN {
    na = split(NUM, a, "\\.")
    nb = split(REF, b, "\\.")
    n = (na < nb ? na : nb)
    for (i = 1 ; i <= n ; i++) {
      if (a[i] < b[i]) exit(1)
      if (a[i] > b[i]) exit(0)
    }
      
    exit(n > 0 ? 0 : 2)
}
