import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.patches import Rectangle
import numpy as np

def make_xaxis(ax, yloc, offset=0.05, **props):
    xmin, xmax = ax.get_xlim()
    locs = [loc for loc in ax.xaxis.get_majorticklocs()
            if loc>=xmin and loc<=xmax]
    tickline, = ax.plot(locs, [yloc]*len(locs),linestyle='',
            marker=lines.TICKDOWN, **props)
    axline, = ax.plot([xmin, xmax], [yloc, yloc], **props)
    tickline.set_clip_on(False)
    axline.set_clip_on(False)
    for loc in locs:
        if loc == 0: loc = 1
        st = ''
        if loc > 500000:
            st = ('%.1f' % (loc / 1000000.)) + 'M'
        elif loc > 10000:
            st = ('%.0f' % (loc / 1000.)) + 'K'
        else:
            st = str(loc)
        st = st.replace('.0', '')
        if loc == 1: loc += 0.005 * (xmax - xmin)
        # if it's close to the edge, bump it a bit.
        ax.text(loc, yloc-offset, '%s'  %st,
                horizontalalignment='center',
                verticalalignment='top')

props = dict(color='black', linewidth=1, markeredgewidth=1)

def genomic_axes(fig, bpmin, bpmax, props=props):
    ax = fig.add_axes((0.0001, 0, .9999, 0.7), frameon=False, yticks=())
    ax.axison = False
    ax.set_autoscale_on(False)

    w = bpmax - bpmin #4100003
    h = 2
    ax.set_xlim(bpmin, bpmax)
    ax.set_ylim(-h, h)

    make_xaxis(ax, 0, offset=0.11, **props)
    r = Rectangle((0, -h/2), w, h, fc=(1, 1, 1, 1))
    ax.add_artist(r)
    return ax

def paint(ax, f, txt):
    if txt is True: txt = f.attribs['ID']
    if txt is False: txt = ''
    strand = f.strand == '-' and -1 or 1
    txt = txt[:11]

    x = (f.start + f.end) / 2
    y = strand * (0.03 if txt else 0.3)

    ax.annotate(txt, xy=(x, y), xycoords='data',
            xytext=(x, strand * 0.94),
                arrowprops=dict(fc='gray', ec=None, width=0.2, headwidth=5),
                #arrowprops=dict(arrowstyle="wedge,tail_width=0.7",
                #                fc="0.6", ec="none"),
                horizontalalignment='right' if f.strand == '-' else 'left',
            verticalalignment='bottom' if f.strand == '-' else 'top',
            fontsize=10,
            rotation=45 if txt else None)



#paint(ax, Feature(2000000, 2000000, '+'), 'at2g25640')
#paint(ax, Feature(3000000, 3000000, '-'), 'bbbat2g52460')
#paint(ax, Feature(1000000, 1000000, '-'), False)

def paint_file(fname, seqid, save_as=None):
    f = plt.figure()

def paint_features(feats, figure, ax=None):
    if ax is None:
        maxbp = max(f.start for f in feats) * 1.1
        ax = genomic_axes(figure, 0, maxbp)

    for f in feats:
        paint(ax,  f, True)


if __name__ == "__main__":
    import sys
    figure = plt.figure()
    from to_gff import to_gff_lines, Feature

    feats = [f for f in to_gff_lines(sys.argv[1], as_obj=True) if f.seqid == '1']
    paint_features(feats, figure)


    plt.show()
    plt.close()
