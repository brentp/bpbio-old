this is a copy of dagchainer:
http://dagchainer.sourceforge.net/


Changes
-------
1. the calling script has been converted to python. dag_chainer.py
   the options should be identical.

2. it removes the need for temporary files (per chr-chr pair) by 
   reading/writing to and from stdout.

3. does less filtering than the perl eqivalent.

4. the Makefile was update to include compile flags that make the resulting
   executable up to 2 times as fast.

5. added a find_nearby.py script. dagchainer finds the *optimal* path through a list
   of hits, but we may also want to include nearby hits that weren't on the optimial
   path. this uses scipy's cKDTree to fast queries to find the close-by points.
   As such scipy>=0.8 is required to use this script.


Example
-------

1. build the dagchainer executable and the python libs ::

  make

2. an example shell script using the atat.dag from arabidopsis ::

    IN=atat.dag
    OUT=a.out
    python dagtools/evalue_adjust.py -c 4 $IN | python dag_chainer.py -g 40000 -D 160000 -A 4 --merge $OUT -
    python dagtools/find_nearby.py -d 20000 --diags ${OUT}.merge --all ${INFILE} > ${OUT}.merge.all
    python dagtools/plot_dag.py -q athaliana_3 -s athaliana_5 -p dotplot.png -d ${OUT}.merge.all


3. and view the output in dotplot.png. the ${OUT}.merge file has the diagonals
   merged, the ${OUT}.merge.all includes points that are near the diagonals, 
   but not on what dagchainer considers the optimal path.
