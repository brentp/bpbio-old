import sys
import blast_misc
import time
import operator

blast_file = 'data/t.blast'
b = blast_misc.blast_array(blast_file, dopickle=0)
print b[:10]
