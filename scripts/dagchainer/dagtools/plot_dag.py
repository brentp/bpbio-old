import sys
from cdagline import DagLine
import numpy as np
from itertools import cycle

from matplotlib import pyplot as plt
f = plt.figure()
ax = f.add_subplot(1, 1, 1)

reverse = False
xmax = 0
ymax = 0
LINES = False

colors = cycle('rgbcmyk')

c = 'y'
pts = []
for line in open(sys.argv[1]):
    if line[0] == '#': 
        if sys.argv[2] in line and sys.argv[3] in line:
            c = colors.next()
        continue
    dag = DagLine(line)

    if dag.a_seqid != sys.argv[2]: continue
    if dag.b_seqid != sys.argv[3]: continue

    if LINES:
        ax.plot([dag.a_start, dag.a_end], 
                [dag.b_start, dag.b_end], c=c)
    else:
        pts.append((dag.a_start, dag.b_start, c))

    if dag.a_end > xmax: xmax = dag.a_end
    if dag.b_end > ymax: ymax = dag.b_end

if not LINES:
    pts = np.array(pts, dtype=[('x', int), ('y', int), ('c', 'S1')])
    ax.scatter(pts['x'], pts['y'], edgecolor='none', c=pts['c'], s=2)


ax.set_xlim(0, xmax)
ax.set_ylim(0, ymax)
plt.show()
plt.close()
