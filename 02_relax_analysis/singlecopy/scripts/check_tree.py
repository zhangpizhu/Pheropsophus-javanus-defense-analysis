import sys
t = open(sys.argv[1]).read()
ok = ("{TEST}" in t) and (";" in t) and (len(t.strip()) > 10)
sys.exit(0 if ok else 1)
