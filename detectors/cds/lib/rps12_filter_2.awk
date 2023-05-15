function min(a,b) {return (a <= b) ? a:b } 
                   (old6 == 1) {
                        print old
                        oldprint = 1
                   }

((old6 == 2 && $6==2) ||
full == 1) {
    print old
    full = 0
}

(((old6 == 2 && $6==3) ||
(old6 == 3 && $6==2)) && full != 1) {
    $1 = old1
    $6 = min(old6,$6)
    full = 1
}

END {print old}

{
    old = $0
    old1 = $1
    old6= $6
}