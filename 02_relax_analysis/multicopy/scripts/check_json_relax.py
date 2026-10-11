import sys, json
# check_json_relax.py v2.7
# v2.7: also require an actual p-value inside "test results" -- a JSON that
# parses but carries no result block must never be resumed as done.
try:
    j = json.load(open(sys.argv[1]))
    tr = j.get("test results", {})
    ok = isinstance(j, dict) and bool(tr) and tr.get("p-value") is not None
    sys.exit(0 if ok else 1)
except Exception:
    sys.exit(1)
