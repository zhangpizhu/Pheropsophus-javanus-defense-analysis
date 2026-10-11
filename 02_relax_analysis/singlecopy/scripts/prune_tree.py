import sys, os

aln_fp, out_fp, og = sys.argv[1], sys.argv[2], sys.argv[3]
tree_str = open(os.environ["SPECIES_TREE"]).read().strip()
test_sp  = os.environ["TEST_SPECIES"]

# species present in this family's alignment
species = set()
for line in open(aln_fp, errors="ignore"):
    if line.startswith(">"):
        species.add(line[1:].split()[0])

# ---- minimal Newick parser (supports bootstrap labels & branch lengths) ----
def parse(s):
    s = s.strip().rstrip(";")
    pos = 0
    def parse_node():
        nonlocal pos
        if s[pos] == "(":
            pos += 1
            children = []
            while True:
                children.append(parse_node())
                if s[pos] == ",":
                    pos += 1
                else:
                    break
            assert s[pos] == ")"
            pos += 1
            label, blen = parse_label_len()
            return {"children": children, "label": label, "len": blen}
        else:
            label, blen = parse_label_len()
            return {"children": None, "label": label, "len": blen}
    def parse_label_len():
        nonlocal pos
        label = ""; blen = ""
        while pos < len(s) and s[pos] not in ",:)":
            label += s[pos]; pos += 1
        if pos < len(s) and s[pos] == ":":
            pos += 1
            while pos < len(s) and s[pos] not in ",()":
                blen += s[pos]; pos += 1
        return label, blen
    return parse_node()

def tips(node):
    if node["children"] is None:
        return [node["label"]]
    out = []
    for c in node["children"]:
        out += tips(c)
    return out

def prune(node, keep):
    # returns (node-or-None, accumulated length)
    if node["children"] is None:
        return (node, 0.0) if node["label"] in keep else (None, 0.0)
    newch = []
    for c in node["children"]:
        cc, _ = prune(c, keep)
        if cc is not None:
            newch.append(cc)
    if not newch:
        return None, 0.0
    if len(newch) == 1:
        ch = newch[0]
        try:
            ch["len"] = str(float(ch["len"] or 0) + float(node["len"] or 0))
        except ValueError:
            pass
        return ch, 0.0
    node["children"] = newch
    return node, 0.0

def serialize(node):
    if node["children"] is None:
        lbl = node["label"]
        if lbl == test_sp:
            lbl = lbl + "{TEST}"
        return lbl + (":" + node["len"] if node["len"] else "")
    inner = ",".join(serialize(c) for c in node["children"])
    return "(" + inner + ")" + (node["label"] or "") + (":" + node["len"] if node["len"] else "")

tree = parse(tree_str)
tree_tips = set(tips(tree))
missing = species - tree_tips
if missing:
    sys.exit(f"[ERROR] {og}: species in alignment but not in tree: {missing}")
keep = tree_tips & species
tree2, _ = prune(tree, keep)
newick = serialize(tree2) + ";"
assert "{TEST}" in newick, f"{og}: {test_sp} not found in this family"
open(out_fp, "w").write(newick + "\n")
print(f"[step1] {og}: pruned to {len(keep)} species, TEST labeled")
