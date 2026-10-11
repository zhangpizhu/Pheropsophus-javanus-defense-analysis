import sys
def read_fa(fp):
    d={}; name=None; buf=[]
    for line in open(fp, errors="ignore"):
        line=line.rstrip("\n\r")
        if line.startswith(">"):
            if name is not None: d[name]="".join(buf)
            name=line[1:].split()[0]; buf=[]
        else:
            buf.append(line.strip())
    if name is not None: d[name]="".join(buf)
    return d
codon=read_fa(sys.argv[1]); aln=read_fa(sys.argv[2])
if not codon: sys.exit(1)
W=max(len(s) for s in aln.values())
ok = all(len(s)==3*W for s in codon.values())
sys.exit(0 if ok else 1)
