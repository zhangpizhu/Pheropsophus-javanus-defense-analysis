import sys
# check_codon_ready.py v2.6
# incremental-reuse gate for 04_pal2nal/<og>.codon.fasta (post-trim mapping).
# v2.6: also reject stop codons (TAA/TAG/TGA) so a stale stop-containing
# codon file is never reused; forces a remap instead.
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
codon=read_fa(sys.argv[1]); trim=read_fa(sys.argv[2])
if not codon: sys.exit(1)
W=max(len(s) for s in trim.values())
ACGT=set("ACGTacgt")
STOPC=("TAA","TAG","TGA")
ok = all(len(s)==3*W for s in codon.values())
usable = sum(1 for s in codon.values() if set(s) & ACGT)
stops = 0
for s in codon.values():
    up=s.upper()
    for j in range(0,len(up)-2,3):
        if up[j:j+3] in STOPC: stops+=1
ok = ok and usable>=4 and usable==len(codon) and stops==0
if not ok:
    sys.stderr.write(f"check_codon_ready: FAIL widths_bad={not all(len(s)==3*W for s in codon.values())} usable={usable}/{len(codon)} stop_codons={stops}\n")
sys.exit(0 if ok else 1)
