import blast_misc
import sys
sys.path.append('/home/bpederse/work/bio/pydiags')
import blastarray
import time
import operator

#blast_file = '/tmp/all.blast'
blast_file = '/tmp/arabidopsischr5_arabidopsischr5.blast'
t0 = time.time()
a = blast_misc.blast_array(blast_file, dopickle=0)
t1 = time.time()
b = blastarray.blast_array(blast_file, dopickle=0)
t2 = time.time()
print t1 - t0
print t2 - t1
print sorted(a[a['query'] == 'AT5G67620'], key=operator.itemgetter('subject'))
print sorted(b[b['query'] == 'AT5G67620'], key=operator.itemgetter('subject'))

print (t2 - t1)/(t1 - t0 * 1.0)
