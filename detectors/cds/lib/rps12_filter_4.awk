function convert1p(p) {
    if (p+0 < L2) {
    I = 1
    if (S1=="F") {
        S = 1
        B = F1
    } else {
        S = -1
        B = T1
    }
    } else {
    I = L2
    if (S2=="F") {
        S = 1
        B = F2
    } else {
        S = -1
        B = T2
    }
    }
    return S*(p - I) + B
}
function convert(p1,p2) {
    p1  = convert1p(p1)
    p2  = convert1p(p2)
    if (p1 < p2)
        res = p1 ".." p2
    else
        res = "complement(" p2 ".." p1 ")"
    return res
}
/[0-9]+\.\.[0-9]+/ {
    s = $0
    r = $0
    while (length(s) > 0) {
        match(s,/[0-9]+\.\.[0-9]+/)
        range = substr(s,RSTART,RLENGTH)
        s = substr(s,RSTART+RLENGTH+1)
        match(range,/^[0-9]+/)
        from = substr(range,RSTART,RLENGTH)
        match(range,/[0-9]+$/)
        to = substr(range,RSTART,RLENGTH)
        sub(range,convert(from,to),r)
    }
    $0=r
    }
{print $0}