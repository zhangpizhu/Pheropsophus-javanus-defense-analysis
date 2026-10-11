import sys
name=None; n=0; cur=0; L=set()
for line in open(sys.argv[1], errors="ignore"):
    line=line.rstrip("\n\r")
    if line.startswith(">"):
        if name is not None: L.add(cur)
        name=line; n+=1; cur=0
    else:
        cur+=len(line.strip())
if name is not None: L.add(cur)
ok = (n>=2) and (len(L)==1) and (next(iter(L),0)>0)
sys.exit(0 if ok else 1)
