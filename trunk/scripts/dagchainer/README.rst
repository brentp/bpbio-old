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


Example
-------

1. build the dagchainer executable and the python libs ::

  make

2. an example shell script using the atat.dag from arabidopsis :

    IN=atat.dag
    OUT=a.out
    python dagtools/evalue_adjust.py -c 4 $IN | python dag_chainer.py -g 40000 -D 160000 -A 4 --merge $OUT -
    python dagtools/plot_dag.py -q athaliana_3 -s athaliana_5 -p dotplot.png -d ${OUT}.merge


3. and view the output in dotplot.png
