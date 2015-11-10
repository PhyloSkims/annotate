#
# NWS alignment
#

# -----------------------
# make identity substitution matrix for given alphabet
#
function IdentityMatrix(alpha, mat, diag, ndiag, _local_, n, i, j, ci, cj) {
  if (diag == "")  diag  = 0
  if (ndiag == "") ndiag = 1
  alpha = toupper(alpha)
  delete mat
  n = length(alpha)
  for (i = 1 ; i <= n ; i++) {
    ci = substr(alpha, i, 1)
    for (j = 1 ; j <= n ; j++) {
      cj = substr(alpha, j, 1)
      mat[ci][cj] = (i == j) ? diag : ndiag
    }
  }
}

# -----------------------
# internal utility : reverse string
#
  function _Reverse(s, _local_, i, n, rs) {
    rs = "";
    n = length(s);
    for (i = n ; i >= 1 ; i--)
      rs = rs "" substr(s, i, 1)
    return rs;
  }

# -----------------------
# internal utility : alignment traceback for NWS
#
  function _Traceback(dir, s1, s2, n1, n2, align, _local_, i, c1, c2, c3) {
    delete align
    while (dir[n1][n2] != 0) {
      if (dir[n1][n2] == "s") {
        c1 = substr(s1, n1--, 1)
        c2 = substr(s2, n2--, 1)
        c3 = (c1 == c2) ? tolower(c1) : "*"
      }
      else if (dir[n1][n2] == "i") {
        c1 = "-"
        c2 = substr(s2, n2--, 1)
        c3 = "-"
      }
      else {
        c1 = substr(s1, n1--, 1)
        c2 = "-"
        c3 = "-"
      }
      align[1] = align[1] "" c1
      align[2] = align[2] "" c2
      align[3] = align[3] "" c3
    }
    for (i = 1 ; i <= 3 ; i++)
      align[i] = _Reverse(align[i])
  }

# -----------------------
# internal utility : min
#
  function _Min(a, b) {
    return (a < b ? a : b)
  }

# -----------------------
# sequence alignment NWS
#
# todo : check alphabet
#
#   --> i
#  |
#  v d
#
function AlignNWS(s1, s2, subst, indel, local, align,
                     _local_, rev, n1, n2, i, j, c1, c2, m2,
                           ws, wi, wd, w,
                           mat, dir) {
                           
   s1 = toupper(s1) ; s2 = toupper(s2)
   n1 = length(s1)  ; n2 = length(s2)
   if (local && (n2 < n1)) {
     rev = s1 ; s1 = s2 ; s2 = rev
     rev = n1 ; n1 = n2 ; n2 = rev
   }

   if (indel == "") indel = 1

   #
   # nws alignment
   #
   for (i = 0 ; i <= n1 ; i++) {
     c1 = substr(s1, i, 1)
     for (j = 0 ; j <= n2 ; j++) {
       c2 = substr(s2, j, 1)
       if (i && j) {
         ws = mat[i-1][j-1] + subst[c1][c2]
         wd = mat[i-1][j] + indel
         wi = mat[i][j-1] + indel
         w  = _Min(ws, _Min(wi, wd))
         mat[i][j] = w
         dir[i][j] = (w == ws) ? "s" : (w == wi) ? "i" : "d"
       } else if (i) {
         mat[i][j] = mat[i-1][j] + indel
         dir[i][j] = "d"
       } else if (j) {
         mat[i][j] = mat[i][j-1] + (local ? 0 : indel)
         dir[i][j] = "i"
       } else {
         mat[i][j] = 0
         dir[i][j] = 0
       }
     }
     # delete mat[i-1]
   }

   #
   # adjust last line in local mode
   #
   if (local) {
     m2 = n2
     for (j = m2 ; j >= 0 ; j--) {
       if (mat[n1][j] < mat[n1][m2])
         m2 = j
     }
     mat[n1][n2] = mat[n1][m2]
     for (j = m2 + 1 ; j <= n2 ; j++) {
       dir[n1][j] = "i"
     }
   }

   #
   # traceback
   #
   _Traceback(dir, s1, s2, n1, n2, align)
   
   if (rev) {
     rev = align[1] ; align[1] = align[2]; align[2] = rev
   }
   
   return mat[n1][n2]
}
