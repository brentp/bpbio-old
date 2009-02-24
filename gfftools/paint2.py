import matplotlib.pyplot as plt
import matplotlib.lines as lines
from matplotlib.patches import Rectangle
from matplotlib.colors import hex2color
import numpy as np

def make_yaxis(ax, xloc, offset=0.05, **props):
    ymin, ymax = ax.get_ylim()
    locs = [loc for loc in ax.yaxis.get_majorticklocs()
            if loc>=ymin and loc<=ymax]
    tickline, = ax.plot([xloc] * len(locs), locs,linestyle='',
            marker=lines.TICKRIGHT, **props)
    axline, = ax.plot([xloc, xloc], [ymin, ymax], **props)
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
        if loc == 1: loc += 0.005 * (ymax - ymin)
        # if it's close to the edge, bump it a bit.
        ax.text(xloc + offset, loc, '%s'  %st,
                fontsize=12,
                zorder=20000,
                color='#333333',
                verticalalignment='center',
                horizontalalignment='left')

props = dict(color='#333333', linewidth=2, markeredgewidth=1)

def genomic_axes(fig, bpmin, bpmax, props=props):
    w = 0.5
    ax = fig.add_axes((0.0001, 0, .9999, .9999), frameon=False, xticks=())
    ax.axison = False
    ax.set_autoscale_on(False)

    h = bpmax - bpmin #4100003
    ax.set_xlim(-w/2, w/2)
    ax.set_ylim(bpmin, bpmax)

    make_yaxis(ax, 0, offset=0.03 * w, **props)
    r = Rectangle((-w/2, 0), w, bpmax, fc=(1, 1, 1, 1))
    ax.add_artist(r)
    return ax

def get_colors(f):
    # let the color be set per feature
    color = f.attribs.get('color', f.attribs.get('fc')) #, 'gray'))
    ec = f.attribs.get('ec', f.attribs.get('edgecolor', f.attribs.get('edge_color')))
    if color:
        if color[0] != '#': color = '#' + color
        color = hex2color(color)
    else: color = 'gray'

    if ec:
        if ec[0] != '#': ec = '#' + ec
        ec = hex2color(ec)
    else:
        ec = color
    return color, ec

def paint(ax, f, txt):
    if txt is True: txt = f.attribs['ID']
    if txt is False: txt = ''
    xmin, xmax = ax.get_xlim()
    strand = f.strand == '-' and -1 or 1
    txt = txt[:11]

    y = (f.start + f.end) / 2
    x = strand * (0.03) #ax.get_xbound()[0]

    color, ec = get_colors(f)

    if f.attribs.get('span','').lower() in ('true', 't', '1', 'yes', 'y'):
        vspan = ax.axhspan(f.start, f.end, xmin=0, xmax=1
                , fill=True
                , facecolor='0.5'
                #, hatch='+'
                , alpha=0.2
                #, fc=(0.5, 0.5, 0.5)
                #, ec=(0,0,0)
                , zorder=10)
        vspan.set_clip_on(True)
    else:
        ax.annotate(txt, xy=(x, y), xycoords='data',
            xytext=(strand * 0.19 * (xmax - xmin), y),
                arrowprops=dict(fc=color, ec=ec, width=1.5, headwidth=5,
                    shrink=0.08),
                #arrowprops=dict(arrowstyle="wedge,tail_width=0.7",
                #                fc="0.6", ec="none"),
                horizontalalignment='right' if f.strand == '-' else 'left',
                verticalalignment='center',
            fontsize=10,
            rotation=None)
            #rotation=45 if txt else None)



#paint(ax, Feature(2000000, 2000000, '+'), 'at2g25640')
#paint(ax, Feature(3000000, 3000000, '-'), 'bbbat2g52460')
#paint(ax, Feature(1000000, 1000000, '-'), False)

def paint_file(fname, seqid, save_as=None):
    f = plt.figure(figsize=(2,10))

def paint_features(feats, figure, ax=None):
    if ax is None:
        maxbp = max(f.start for f in feats) * 1.1
        ax = genomic_axes(figure, 0, maxbp)

    for f in feats:
        paint(ax,  f, True)

    

if __name__ == "__main__":
    import sys
    figure = plt.figure(figsize=(2, 10))
    from to_gff import to_gff_lines, Feature

    feats = [f for f in to_gff_lines(sys.argv[1], as_obj=True) if f.seqid == sys.argv[2]]
    paint_features(feats, figure)


    plt.show()
    plt.close()
