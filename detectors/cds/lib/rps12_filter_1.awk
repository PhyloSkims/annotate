function overlap(x1,y1,x2,y2) {  
                      return (((x1+0 <= x2+0) && ((y1+1) >= x2+0)) || 
                              ((x2+0 <= x1+0) && ((y2+1) >= x1+0))) 
                    } 
function min(a,b) {return (a <= b) ? a:b } 
function max(a,b) {return (a >= b) ? a:b }

(NR==1) {i=0
            frg[i]=$0 
        } 

(x1 && y1) { 
    if (overlap(x1,y1,$1,$2)) {
        $1 = min(x1,$1) 
        $2 = max(y1,$2) 
        if (overlap(v1,w1,$3,$4)) { 
            $3 = min(v1,$3) 
            $4 = max(w1,$4) 
            } 
        } 
        else i++ 
} 

(x1 && y1) { 
    frg[i] = $0
} 

{ x1 = $1 
    y1 = $2 
    v1 = $3 
    w1 = $4 
} 

END {
    for (j = 0; j <= i; j++) {
        print frg[j]
    }
}