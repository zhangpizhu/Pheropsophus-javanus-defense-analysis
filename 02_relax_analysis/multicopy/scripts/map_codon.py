import sys
# map_codon.py v2.4/v2.5: trim-set keyed + drop no-ACGT sequences
aa_raw_fp, aa_trim_fp, codon_fp, out_fp = sys.argv[1:5]

def read_fasta(fp):
    order = []; d = {}; name = None; buf = []
    for line in open(fp, errors="ignore"):
        line = line.rstrip("\n\r")
        if line.startswith(">"):
            if name is not None: d[name] = "".join(buf)
            name = line[1:].split()[0]; order.append(name); buf = []
        else:
            buf.append(line.strip())
    if name is not None: d[name] = "".join(buf)
    return order, d

order_raw, aa_raw = read_fasta(aa_raw_fp)
order_trim, aa_trim = read_fasta(aa_trim_fp)
_, codon = read_fasta(codon_fp)
missing_trim = [k for k in order_trim if k not in aa_raw]
missing_cod  = [k for k in order_trim if k not in codon]
assert not missing_trim, f"trim sequences not in raw alignment: {missing_trim[:5]}"
assert not missing_cod,  f"trim sequences not in pal2nal output: {missing_cod[:5]}"
W = len(aa_raw[order_raw[0]])
assert all(len(aa_raw[k]) == W for k in order_raw), "AA alignment widths differ"
assert all(len(codon[k]) == 3*W for k in order_raw), "codon length != 3 x AA width"
use_order = order_trim
rrows = [aa_raw[k] for k in use_order]
trows = [aa_trim[k] for k in use_order]
n_trim = len(trows[0])
assert all(len(s) == n_trim for s in trows), "trimmed alignment widths differ"
raw_sig = ["".join(r[i] for r in rrows) for i in range(W)]
trim_sig = ["".join(t[j] for t in trows) for j in range(n_trim)]
kept = []
j = 0
for i in range(W):
    if j < n_trim and raw_sig[i] == trim_sig[j]:
        kept.append(i); j += 1
assert j == n_trim, f"trimAL column matching failed at {j}/{n_trim}"
assert len(kept) > 0, "no columns retained"
ACGT = set("ACGTacgt")
dropped_nodata = []
with open(out_fp, "w") as o:
    for k in use_order:
        seq = codon[k]
        if not (set(seq) & ACGT):
            dropped_nodata.append(k); continue
        o.write(f">{k}\n" + "".join(seq[3*c:3*c+3] for c in kept) + "\n")
print(f"[step4] {out_fp}: kept {len(kept)} codon columns / {W} cols; "
      f"{len(order_raw)-len(order_trim)} dropped-by-trimal; "
      f"{len(dropped_nodata)} no-ACGT dropped")
assert len(use_order) - len(dropped_nodata) >= 4, "fewer than 4 usable sequences"
