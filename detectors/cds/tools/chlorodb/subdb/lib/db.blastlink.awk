#

function min(x, y) {
  return ((x < y) ? x : y)
}

BEGIN {
  if (COVMIN == "") COVMIN = 50
  if (PMAX == "")   PMAX   = 1e-6
  if (IDMIN == "")  IDMIN  = 30
}

/^#/ {
  hitnum = 0;
  next;
}

{
  if ($1 == $2) next
  
  hitnum++;
  
  na = split($1, a, "@");
  if (na < 2) {
    print "query file not properly formatted" > "/dev/stderr"
    exit(1);
  }
  len1  = a[na];

  na = split($2, a, "@");
  if (na < 2) {
    print "bank file not properly formatted" > "/dev/stderr"
    exit(1);
  }
  len2 = a[na];
  
  id  = $3 + 0.0;
  ali = $4;

  covmin = ali * 100. / min(len1, len2);
  
  proba = $11 + 0.0;
  
  if ((covmin > COVMIN) && ((proba < PMAX) || (proba == 0)) && (id > IDMIN)) {
    print $1, $2, hitnum, id, covmin, proba, ali, len1, len2;
  }
}

