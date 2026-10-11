import sys, os

# v2.4 helper: replace stop codons (TAA/TAG/TGA, case-insensitive) with '---'
# so HyPhy RELAX accepts the input alignment.
#
# Root cause (110 single-copy run, 2026-09-09):
#   HyPhy RELAX asserts the input alignment must not contain stop codons
#   ("assert(relax.codon_filter.sites*3==relax.codon_data.sites)").
#   26/42 control-B families each aborted at startup because ONE sequence in
#   the alignment carried stop codon(s), mostly a retained terminal stop.
#   Scattered species: ILU x4, PCH x4, LSP x3, PPY x3, APA x2, NSA x2,
#   CFO x2, OAR x2, PNI x2, CPR/AFU/PCH... one bad sequence is enough to
#   kill the whole family.
# Why '---' (gaps) and not deletion:
#   Per-sequence codon deletion would desync the alignment columns; masking
#   keeps the alignment rectangular and treats those codons as missing data.
#   Gaps are proven-safe for this HyPhy install (the 106 multicopy pipeline
#   ran gapped codon alignments -> 238 finished RELAX jsons).
#   All-gap columns (only possible if masking empties a whole column) are
#   dropped afterwards - they carry zero information and keep the HyPhy
#   codon filter happy.
# usage: python3 fix_stops.py IN.fasta OUT.fasta
# exit 0 = sanitized (even if zero stops found); exit != 0 = data problem.

STOPS = {"TAA", "TAG", "TGA"}

def read_fasta(fp):
    names, seqs = [], []
    cur, buf = None, []
    for line in open(fp, errors="ignore"):
        line = line.rstrip("\r\n")
        if line.startswith(">"):
            if cur is not None:
                seqs.append("".join(buf))
            cur = line[1:].rstrip()
            names.append(cur)
            buf = []
        elif line.strip():
            buf.append(line.strip())
    if cur is not None:
        seqs.append("".join(buf))
    return names, seqs

def main():
    in_fp, out_fp = sys.argv[1], sys.argv[2]
    names, seqs = read_fasta(in_fp)
    if not names:
        sys.exit("[ERROR] no sequences found in " + in_fp)
    L = len(seqs[0])
    if L == 0:
        sys.exit("[ERROR] first sequence is empty in " + in_fp)
    if L % 3 != 0:
        sys.exit("[ERROR] alignment length %d not divisible by 3" % L)
    total_stops = 0
    masked_seqs = []
    for nm, s in zip(names, seqs):
        tag = nm.split()[0] if nm.split() else "?"
        if len(s) != L:
            sys.exit("[ERROR] %s: length %d != %d (alignment not rectangular)" % (tag, len(s), L))
        chunks = []
        n_stop = 0
        for i in range(0, L, 3):
            cod = s[i:i + 3]
            if cod.upper() in STOPS:
                chunks.append("---")
                n_stop += 1
            else:
                chunks.append(cod)
        t = "".join(chunks)
        if n_stop:
            print("  masked %d stop codon(s) in '%s'" % (n_stop, tag))
        total_stops += n_stop
        masked_seqs.append(t)
    # drop all-gap columns (zero information, avoids filter edge cases)
    keep_cols = [i for i in range(L) if any(s[i] not in "-" for s in masked_seqs)]
    n_drop = L - len(keep_cols)
    if n_drop:
        print("  dropped %d all-gap column(s)" % n_drop)
        masked_seqs = ["".join(s[i] for i in keep_cols) for s in masked_seqs]
        L2 = len(keep_cols)
    else:
        L2 = L
    # warn on nearly-empty sequences (data-layer smell, do not auto-drop)
    for nm, s in zip(names, masked_seqs):
        amb = sum(1 for c in s if c == "-" or c.upper() == "N")
        if amb / L2 > 0.5:
            tag = nm.split()[0] if nm.split() else "?"
            print("  WARNING: '%s' is %.0f%% gaps/ambiguous after masking" % (tag, 100.0 * amb / L2))
    tmp = out_fp + ".tmp"
    with open(tmp, "w") as o:
        for nm, s in zip(names, masked_seqs):
            o.write(">" + nm + "\n")
            for i in range(0, len(s), 60):
                o.write(s[i:i + 60] + "\n")
    os.replace(tmp, out_fp)
    print("[fix_stops] %d seqs x %d nt | stop codons masked: %d | all-gap columns dropped: %d" % (len(names), L2, total_stops, n_drop))
    print("[fix_stops] wrote %s" % out_fp)

main()
