"""
take any number of non-sorted gff3 files and spit out a sorted
version that genometools will accept.
my need to add more top-level features to sortable_type()
"""
import sys
import re

files = sys.argv[1:]

lines = []

def sortable_type(ft):
    if ft == 'gene' or ft == 'pseudogene' or ft == 'MIR': return 0
    if ft == 'mRNA': return 1
    if ft == 'exon':  return 3
    if ft == 'CDS':  return 4
    return 2

id_re = re.compile("ID=([^;]+)")
def to_id(l):
    try:
        return re.search(id_re, l).groups(0)[0]
    except:
        #for lines that dont have ID
        # return a unique thing.
        return to_id.counter
    finally:
        to_id.counter += 1
to_id.counter = 0

for fi in files:
    if fi == "-": fi = sys.stdin
    else: fi = open(fi)
    for line in fi:
        if line[0] == '#': continue
        sline = line.split("\t")
        lines.append((sline[0], sortable_type(sline[2]), int(sline[3]), line))

lines.sort()
print '##gff-version 3'
last = lines[-1]
for line in lines:
    # since theyre sorted, can remove dups without an extra 
    # dict.
    sline = line[-1].split()
    if sline == last: continue

    # repeated start, stop, type, ID with something
    # else different.
    if sline[2:5] == last[2:5]:
        if to_id(sline[-1]) == to_id(last[-1]):
            # could also just give it a different id.
            continue

    print line[-1],
    last = sline

"""
from __future__ import with_statement
from tempfile import NamedTemporaryFile
import gt
with NamedTemporaryFile() as tf:
    tf.write('##gff-version 3\n')
    for line in lines:
        tf.write(line[-1])
    tf.flush()
    ins = gt.GFF3InStream(tf.name)
    outs = gt.GFF3OutStream(ins)
    feature = outs.next_tree()
    while feature:
        feature = outs.next_tree()
"""
