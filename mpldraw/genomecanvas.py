from matplotlib.patches import FancyArrow, Rectangle
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from cStringIO import StringIO
from matplotlib.ticker import Formatter, LinearLocator
from mpl_toolkits.axes_grid import axes_size as Size
from mpl_toolkits.axes_grid import Divider

BOTTOM_AXIS, MIDDLE_AXIS, TOP_AXIS = range(3)

class BasePairTickFormatter(Formatter):
    def __call__(self, x, pos=None):
        fmts = ((1e9, 'G'), (1e6, 'M'), (1e3, 'K'))
        for div, fmt in fmts:
            if x * 2 > div:
                x = "%.3f" % (float(x) / div)
                return x.rstrip('0').rstrip('.') + fmt
        return str(int(x))

class BasePairLocator(LinearLocator):
    pass

 

def xy_dxy(start, stop, strand):
        if strand in ('-1', '-', -1):
            x = stop
            dx = start - stop
        else:
            x = start
            dx = stop - start
        # set y to 0.25 as the lables go above.
        return x, 0.25, dx

class Gene(FancyArrow):
    def __init__(self, text, start, stop, strand=None, width=0.14, ec='none', **kwargs):
        assert start <= stop
        self.text = text

        x, y, dx = xy_dxy(start, stop, strand)
        head_length = self.get_head_length(start, stop, strand, kwargs)

        self._textx = start
        self._texty = y + (1.01 * width) + .05
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
    def __init__(self, text, start, stop, strand=None, width=0.25, ec='none', y=None, **kwargs):
        x, _, dx = xy_dxy(start, stop, strand)
        if y is None:
            y = 0.25
        self._textx = start + abs(dx)/ 20.
        self._texty = y
        self.text = text
        w = width
        if not ec in ('none', None) and not kwargs.get('fc'):
            kwargs['fc'] = ec
        Rectangle.__init__(self, (x, y - w / 2.0), width=dx, height=w, ec=ec, **kwargs)

class HSP(Block):
    def __init__(self, start, stop, strand=None, width=0.18, fc='#cacaca', ec='#777777', **kwargs):
        if not 'y' in kwargs: kwargs['y'] = 0.5
        Block.__init__(self, None, start, stop, strand, width=width, fc=fc, ec=ec, **kwargs)

    

class CDS(Block):
    def __init__(self, start, stop, strand=None, **kwargs):
        fc = kwargs.get('fc', 'red')
        kwargs['fc'] = fc
        kwargs['width'] = kwargs.get('width', 0.18)
        kwargs['ec'] = 'black'

        Block.__init__(self, None, start, stop, strand, **kwargs)


class ThreeAxFigure(Figure):
    def __init__(self, figsize=(768, 196), dpi=96):
        dpi = float(dpi)
        figsize = (figsize[0] / dpi, figsize[1] / dpi)
        Figure.__init__(self, figsize=figsize, dpi=dpi)
        self.patch.set_alpha(0.0)

        self.canvas = FigureCanvas(self)
        self.setup_axes()

    def setup_axes(self):

        rect = (0, 0.07, 1, 1) 
        axes = [self.add_axes(rect, autoscale_on=False, aspect='auto', alpha=1.0, label=str(i)) for i in range(3)]

        hori = [Size.AxesX(axes[0])]
        vert = [Size.Scaled(0.20), Size.Scaled(0.60), Size.Scaled(0.20)]

        d = Divider(self, rect, hori, vert) 

        for i, ax in enumerate(axes):
            ax.set_axes_locator(d.new_locator(nx=0, ny=i))

        self._setup_axes()


    def _setup_axes(self):
        """ set the ranges and spines"""

        for ax in self.axes:
            for loc, sp in ax.spines.iteritems():
                if loc != 'bottom': sp.set_color('none')
        self.axes[BOTTOM_AXIS].xaxis.set_major_formatter(BasePairTickFormatter())

        for ax in self.axes:
            ax.set_yticks([])
            ax.patch.set_alpha(0.0)

        for ax in (MIDDLE_AXIS, TOP_AXIS):
            self.axes[ax].set_xticks([])

        for ax in self.axes:
            ax.set_ylim(0, 1)

    def set_xlim(self, xmin, xmax):
        for ax in self.axes:
            ax.set_xlim(xmin, xmax)


    def add_patch(self, patch, axis=None):
        self.axes[axis].add_patch(patch)
        if hasattr(patch, 'text'):
            if not patch.text is None:
                self.axes[axis].text(patch._textx, patch._texty, patch.text, va='center', ha='left', fontsize=8.5)

    def save(self, filename=None, fontsize=8, autoscale=True):

        for xt in self.axes[BOTTOM_AXIS].get_xticklabels():
            xt.set_fontsize(fontsize)
        if autoscale:
            self.axes[MIDDLE_AXIS].axes.autoscale_view(scalex=False, scaley=True)

        if filename is None:
            s = StringIO()
            self.canvas.print_figure(s, dpi=self.dpi)
            s.seek(0)
            return s.getvalue()

        else:
            self.canvas.print_figure(filename, dpi=self.dpi)
    savefig = save


class GenomeFigure(ThreeAxFigure):
    pass


if __name__ == "__main__":

    gf = GenomeFigure()

    g1 = Gene("G1", 0, 100, 1)
    gf.add_patch(g1, TOP_AXIS)

    g2 = Gene("G@", 10, 130, -1)
    gf.add_patch(g2, BOTTOM_AXIS)

    cds = CDS(60, 90, -1)
    gf.add_patch(cds, BOTTOM_AXIS)


    g3 = Gene("At2g25640", 140, 290, -1)
    gf.add_patch(g3, BOTTOM_AXIS)

    g4 = Gene("At2g25640", 1400, 1590, 1)
    gf.add_patch(g4, TOP_AXIS)

    nostrand = Block("someblock", 1000, 1200, fc='red', ec='black', y=0.7)
    gf.add_patch(nostrand, MIDDLE_AXIS)

    gf.add_patch(HSP(1460, 1503), MIDDLE_AXIS)
    import numpy as np
    gf.axes[1].plot(np.sin(np.linspace(0, 10, 1600)))

    gf.set_xlim(0, 1600)
    gf.save('/var/www/t/t.png')

 
    # the buffer
    print len(gf.save())
    #gf.clear()

