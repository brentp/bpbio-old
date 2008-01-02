import blast_misc
import sys
import time
import operator

#blast_file = '/tmp/all.blast'
blast_file = '/tmp/arabidopsischr5_arabidopsischr5.blast'
t0 = time.time()
b = blast_misc.blast_array(blast_file, dopickle=0)
t1 = time.time()
print t1 - t0
print sorted(b[b['query'] == 'AT5G67620'], key=operator.itemgetter('subject'))
