import sys, json
try:
    j = json.load(open(sys.argv[1]))
    sys.exit(0 if isinstance(j, dict) and len(j) > 0 else 1)
except Exception:
    sys.exit(1)
