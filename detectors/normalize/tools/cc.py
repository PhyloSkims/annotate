#!/usr/bin/env python 

import sys

data = open(sys.argv[1])

ccs = []
for line in data:
    
    parts = line.strip().split()
    if len(parts) >= 2:
        a = parts[0]
        b = parts[1]
    else:
        continue
    
    newcc=set([a,b])
    
    keep=set()
    found=set()
    for i in range(len(ccs)):
        if len(found) < 2:
            cc=ccs[i]
            if a not in found and a in cc:
                found.add(a)
                keep.add(i)
            if b not in found and b in cc:
                found.add(b)
                keep.add(i)
                
    for i in keep:
        newcc |= ccs[i]
        
    newccs=[newcc]
                    
    for i in range(len(ccs)):
        if i not in keep:
            newccs.append(ccs[i])
            
    ccs=newccs
    
ccs.sort(key=len, reverse=True)

for i in range(len(ccs)):
    cc=ccs[i]
    for l in cc:
        sys.stdout.write("%d %s\n" % (i,l))
    