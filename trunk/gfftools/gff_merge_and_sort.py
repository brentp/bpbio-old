"""
take any number of non-sorted gff3 files and spit out a sorted
version that genometools will accept.
my need to add more top-level features to sortable_type()
"""
import sys

files = sys.argv[1:]

lines = []

def sortable_type(ft):
    if ft == 'gene' or ft == 'pseudogene': return 0
    if ft == 'mRNA': return 1
    if ft == 'CDS':  return 2
    return 3


for fi in files:
    for line in open(fi):
        if line[0] == '#': continue
        sline = line.split("\t")
        lines.append((sline[0], sortable_type(sline[2]), int(sline[3]), line))

lines.sort()
print '##gff-version 3'
for line in lines:
    print line[-1],

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
