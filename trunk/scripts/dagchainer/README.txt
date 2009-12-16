this is a copy of dagchainer:
http://dagchainer.sourceforge.net/

but the calling script has been converted to python. dag_chainer.py
the options should be identical.
it removes the need for 1 temporary file (per chr-chr pair) by reading
from stdout.

the Makefile was update to include compile flags that make the resulting
executable up to 2 times as fast.
