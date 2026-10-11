import os, json, glob

workdir = os.environ.get("WORKDIR", "RELAX_singlecopy_42")
rows = []
for fp in sorted(glob.glob(os.path.join(workdir, "02_relax_json", "*.RELAX.json"))):
    og = os.path.basename(fp).replace(".RELAX.json", "")
    try:
        j = json.load(open(fp))
        tr = j.get("test results", {})
        k = tr.get("relaxation or intensification parameter", "")
        p = tr.get("p-value", "")
        lr = tr.get("LR", "")
        rows.append((og, k, lr, p))
    except Exception as e:
        rows.append((og, "PARSE_ERROR", "", str(e)[:60]))
out = os.path.join(workdir, "RELAX_summary.tsv")
with open(out, "w") as o:
    o.write("OG\tk\tLR\tp_value\n")
    for r in rows:
        o.write("\t".join(str(x) for x in r) + "\n")
print(f"[summary] {len(rows)} families -> {out}")
