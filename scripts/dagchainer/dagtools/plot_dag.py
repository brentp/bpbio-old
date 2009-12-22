import sys
from cdagline import DagLine
import numpy as np
from itertools import cycle
import matplotlib

def main(args):

    import optparse
    p = optparse.OptionParser("plot a dagfile")
    p.add_option('-q', '--qseqid', dest='qseqid', help="seqid of the query")
    p.add_option('-s', '--sseqid', dest='sseqid', help="seqid of the subject")
    p.add_option('-p', '--png', dest='png', help="path of png to save file.")

    p.add_option('-l', '--lines', dest='lines', help="draw as lines (not dots)",
                 action='store_true')

    p.add_option('-d', '--dag', dest='dag', help='path to dag file')
    # TODO outfile.

    opts, _ = p.parse_args(args)
    if not (opts.qseqid and opts.sseqid and opts.dag and opts.png):
        sys.exit(p.print_help())
    plot(opts.dag, opts.qseqid, opts.sseqid, opts.png, opts.lines)

def plot(dagfile, qseqid, sseqid, png, lines=False):

    matplotlib.use('Agg')

    from matplotlib import pyplot as plt
    f = plt.figure()
    ax = f.add_subplot(1, 1, 1)

    xmax = 0
    ymax = 0

    colors = cycle('rgbcmyk')

    c = 'y'
    pts = []
    for line in open(dagfile):
        if line[0] == '#': 
            if qseqid in line and sseqid in line:
                c = colors.next()
            continue
        dag = DagLine(line)

        if dag.a_seqid != qseqid: continue
        if dag.b_seqid != sseqid: continue

        if lines:
            ax.plot([dag.a_start, dag.a_end], 
                    [dag.b_start, dag.b_end], c=c)
        else:
            pts.append((dag.a_start, dag.b_start, c))

        if dag.a_end > xmax: xmax = dag.a_end
        if dag.b_end > ymax: ymax = dag.b_end

    if not lines:
        pts = np.array(pts, dtype=[('x', int), ('y', int), ('c', 'S1')])
        ax.scatter(pts['x'], pts['y'], edgecolor='none', c=pts['c'], s=2)


    ax.set_xlim(0, xmax)
    ax.set_ylim(0, ymax)
    plt.savefig(png)
    #plt.show()
    plt.close()

if __name__ == "__main__":
    if len(sys.argv) < 2: sys.argv.append('zzz')
    main(sys.argv[1:])
