import sys
best = {}
for line in open(sys.argv[1]):
    line = line[:-1].split("\t")
    key = (line[0], line[1])
    if key in best:
        best[key].append([float(line[-1])])
    else: best[key] = [float(line[-1])]

b2 = {}
for k,v in best.iteritems():
    v.sort()
    b2[k] = v[:(len(v) < 4 and len(v) or 3)]

for li in open(sys.argv[1]):
    line = li[:-1].split("\t")
    key = (line[0], line[1])
    if float(line[-1]) in b2[key]:
        print li,


