from matplotlib.patches import FancyArrow
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from cStringIO import StringIO

class Gene(FancyArrow):
    def __init__(self, text, start, stop, strand, width=0.07, ec='none', **kwargs):
        assert start < stop
        self.text = text

        dy = 0
        if strand in ('-1', '-', -1):
            x = stop
            dx = start - stop
            y = -0.22
        else:
            x = start
            dx = stop - start
            y = 0.22

        # move it down a bit since text goes above.
        y -= 0.05

        head_length = (stop - start)/ 4
        #if head_length > 200: head_length = 200
        if head_length < 10: head_length = max(10, stop - start)
        self._textx = start
        self._texty = y + (1.8 * width) + .16
        FancyArrow.__init__(self, x, y, dx, dy, head_length=head_length, width=width, length_includes_head=True, ec='none', **kwargs)




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

    def set_xlim(self, xmin, xmax):
        self.ax.set_xlim(xmin, xmax)




if __name__ == "__main__":

    gf = GeneFigure()

    g1 = Gene("G1", 0, 100, 1)
    gf.add_patch(g1)

    g2 = Gene("G@", 10, 130, -1)
    gf.add_patch(g2)


    g3 = Gene("At2g25640", 140, 290, -1)
    gf.add_patch(g3)

    g4 = Gene("At2g25640", 1400, 1590, 1)
    gf.add_patch(g4)

    gf.set_xlim(0, 1600)
    gf.save('t.png')
 
    # the buffer
    print len(gf.save())
