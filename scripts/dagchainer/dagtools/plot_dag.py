import sys
from cdagline import DagLine
from itertools import cycle

from matplotlib import pyplot as plt
f = plt.figure()
ax = f.add_subplot(1, 1, 1)

reverse = False
xmax = 0
ymax = 0

colors = cycle('rgbcky')


for line in open(sys.argv[1]):
    if line[0] == '#': 
        c = colors.next()
        reverse = 'reverse' in line
        continue
    dag = DagLine(line)

    if dag.a_seqid != sys.argv[2]: continue
    if dag.b_seqid != sys.argv[3]: continue

    ax.plot([dag.a_start, dag.a_end], 
            [dag.b_start, dag.b_end], c=c)
    if dag.a_end > xmax: xmax = dag.a_end
    if dag.b_end > ymax: ymax = dag.b_end

ax.set_xlim(0, xmax)
ax.set_ylim(0, ymax)
plt.show()
plt.close()
