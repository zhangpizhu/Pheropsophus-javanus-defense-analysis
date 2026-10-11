import sys
# check_codon_final.py v2.6
# v2.6: reject alignments containing stop codons (TAA/TAG/TGA).
# hyphy relax hard-fails on stop codons; catch it here so the family
# fails with a clear 'codon_invalid' marker instead of a hyphy error.
name=None; n=0; usable=0; cur=''; L=set(); stops=0
STOPC=("TAA","TAG","TGA")
def flush(cur):
    ln=len(cur)
    u=1 if (ln>0 and (set(cur)&set("ACGTacgt"))) else 0
    st=0
    if ln%3==0:
        up=cur.upper()
        for j in range(0,ln,3):
            if up[j:j+3] in STOPC: st+=1
    return ln,u,st
for line in open(sys.argv[1], errors="ignore"):
    line=line.rstrip("\n\r")
    if line.startswith(">"):
        if name is not None:
            ln,u,st=flush(cur); L.add(ln); usable+=u; stops+=st
        name=line; n+=1; cur=''
    else:
        cur+=line.strip()
if name is not None:
    ln,u,st=flush(cur); L.add(ln); usable+=u; stops+=st
ok = (n>=4) and (usable==n) and (len(L)==1) and (next(iter(L),0)>0) and (next(iter(L),0)%3==0) and (stops==0)
if not ok:
    sys.stderr.write(f"check_codon_final: FAIL n={n} usable={usable} widths={sorted(L)} stop_codons={stops}\n")
sys.exit(0 if ok else 1)
