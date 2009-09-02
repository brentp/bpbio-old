from matplotlib.patches import FancyArrow, Rectangle
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from cStringIO import StringIO
from matplotlib.ticker import Formatter, LinearLocator

class BasePairTickFormatter(Formatter):
    def __call__(self, x, pos=None):
        fmts = ((1e9, 'G'), (1e6, 'M'), (1e3, 'K'))
        for div, fmt in fmts:
            if x * 2 > div:
                x = "%.3f" % (float(x) / div)
                x = x.rstrip('0').rstrip('.') + fmt
                break
        else:
            x = str(int(x))
        return x

class BasePairLocator(LinearLocator):
    pass
    

def xy_dxy(start, stop, strand):
        if strand in ('-1', '-', -1):
            x = stop
            dx = start - stop
            y = -0.22
        else:
            x = start
            dx = stop - start
            y = 0.22 if strand else 0.0
        return x, y, dx

class Gene(FancyArrow):
    def __init__(self, text, start, stop, strand=None, width=0.07, ec='none', **kwargs):
        assert start <= stop
        self.text = text

        x, y, dx = xy_dxy(start, stop, strand)
        # move it down a bit since text goes above.
        y -= 0.05
        head_length = self.get_head_length(start, stop, strand, kwargs)

        self._textx = start
        self._texty = y + (0.5 * width) + .16
        if not ec in ('none', None) and not kwargs.get('fc'):
            kwargs['fc'] = ec
        FancyArrow.__init__(self, x, y, dx, 0, head_width=width * 1.35, head_length=head_length, width=width, length_includes_head=True, ec=ec, **kwargs)

    def get_head_length(self, start, stop, strand, kwargs):
        if 'head_length' in kwargs:
            return float(kwargs['head_length'])
        if not strand:
            return 0

        head_length = (stop - start)/ 6.
        #if head_length > 200: head_length = 200
        if head_length < 10: head_length = max(10, stop - start)
        return head_length

class Block(Rectangle):
    def __init__(self, text, start, stop, strand=None, width=0.3, ec='none', **kwargs):
        x, y, dx = xy_dxy(start, stop, strand)
        self._textx = start + abs(dx)/ 20.
        self._texty = y 
        self.text = text
        w = width
        y -= 0.05
        if not ec in ('none', None) and not kwargs.get('fc'):
            kwargs['fc'] = ec
        Rectangle.__init__(self, (x, y - w / 2.0), width=dx, height=w, ec=ec, **kwargs)

class CDS(Block):
    def __init__(self, start, stop, strand=None, **kwargs):
        fc = kwargs.get('fc', 'red')
        kwargs['fc'] = fc
        kwargs['width'] = kwargs.get('width', 0.09)
        kwargs['ec'] = 'black'

        Block.__init__(self, None, start, stop, strand, **kwargs)


class SimpleFigure(Figure):
    def __init__(self, figsize=(512, 256), dpi=96):
        dpi = float(dpi)
        figsize = (figsize[0] / dpi, figsize[1] / dpi)
        Figure.__init__(self, figsize=figsize, dpi=dpi)
        self.patch.set_alpha(0.0)

        self.canvas = FigureCanvas(self)
        ax = self.add_axes([0, 0.22, 1.0, .78], autoscale_on=False, aspect='auto', alpha=0.0)

        for loc, sp in ax.spines.iteritems():
            if loc != 'bottom': sp.set_color('none')
            #else: sp.set_position('center')

        ax.set_yticks([])
        ax.patch.set_alpha(0.0)
        self.ax = ax


    def set_bbox(self, bbox):
        self.ax.set_xlim(bbox[0], bbox[2])
        self.ax.set_ylim(bbox[1], bbox[3])


    def add_patch(self, patch):
        self.ax.add_patch(patch)
        if hasattr(patch, 'text'):
            if not patch.text is None:
                self.ax.text(patch._textx, patch._texty, patch.text, va='top', ha='left', fontsize=8.5)

    def save(self, filename=None):

        for xt in self.ax.get_xticklabels():
            xt.set_fontsize(8)
        if filename is None:
            s = StringIO()
            self.canvas.print_figure(s, dpi=self.dpi)
            s.seek(0)
            return s.getvalue()

        else:
            self.canvas.print_figure(filename, dpi=self.dpi)
    savefig = save


class GeneFigure(SimpleFigure):
    def __init__(self, figsize=(768, 96), dpi=96):
        SimpleFigure.__init__(self, figsize, dpi)
        self.ax.set_ylim(-0.5, 0.5)
        self.ax.xaxis.set_major_formatter(BasePairTickFormatter())
        #self.ax.xaxis.set_major_locator(BasePairLocator())
        #self.ax.grid()

    def set_xlim(self, xmin, xmax):
        self.ax.set_xlim(xmin, xmax)




if __name__ == "__main__":

    gf = GeneFigure()

    g1 = Gene("G1", 0, 100, 1)
    gf.add_patch(g1)

    g2 = Gene("G@", 10, 130, -1)
    gf.add_patch(g2)

    cds = CDS(60, 90, -1)
    gf.add_patch(cds)


    g3 = Gene("At2g25640", 140, 290, -1)
    gf.add_patch(g3)

    g4 = Gene("At2g25640", 1400, 1590, 1)
    gf.add_patch(g4)

    nostrand = Block("someblock", 1000, 1200, fc='red', ec='black')
    gf.add_patch(nostrand)

    gf.set_xlim(0, 1600)
    gf.save('t.png')
 
    # the buffer
    print len(gf.save())
