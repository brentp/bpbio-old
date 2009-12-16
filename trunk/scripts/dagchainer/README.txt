this is a copy of dagchainer:
http://dagchainer.sourceforge.net/


Changes
-------
1. the calling script has been converted to python. dag_chainer.py
   the options should be identical.
2. it removes the need for 1 temporary files (per chr-chr pair) by 
   reading/writing to and from stdout.
3. does less filtering than the perl eqivalent.

4. the Makefile was update to include compile flags that make the resulting
   executable up to 2 times as fast.
