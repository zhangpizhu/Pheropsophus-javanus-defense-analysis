import sys, re
src, dst, prefix = sys.argv[1], sys.argv[2], sys.argv[3]
tree = open(src).read().strip()
labeled = re.sub(r'(' + re.escape(prefix) + r'[^,():\s]*)(?=:)', r'\1{TEST}', tree)
n = labeled.count('{TEST}')
if n == 0:
    tip_names = set(re.findall(r'[(,]([^,():\s]+):', tree))
    has_pja = any(t.startswith(prefix.rstrip('|')) for t in tip_names)
    if has_pja:
        sys.exit(f"TEST labeling failed: {prefix}* tips exist but regex did not match")
    else:
        sys.exit(f"no {prefix} tips in this tree (PJA sequence(s) dropped by codon-mapping no-ACGT filter, or absent in family)")
open(dst, "w").write(labeled + "\n")
print(f"[step6] {src}: labeled {n} TEST branches")
