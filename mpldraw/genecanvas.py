from matplotlib.patches import FancyArrow
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure

class Gene(FancyArrow):
    def __init__(self, text, start, stop, strand, width=0.05, **kwargs):
        assert start < stop
        self.text = text

        dy = 0
        if strand in ('-1', '-', -1):
            x = stop
            dx = start - stop
            y = -0.25
        else:
            x = start
            dx = stop - start
            y = 0.25

        head_length = min(10, max(20, (stop - start)/ 4))
        self._textx = start
        self._texty = y + (1.8 * width) + .16
        FancyArrow.__init__(self, x, y, dx, dy, head_length=head_length, width=width, length_includes_head=True, **kwargs)




class SimpleFigure(Figure):
    def __init__(self, figsize=(512, 256), dpi=96):
        dpi = float(dpi)
        figsize = (figsize[0] / dpi, figsize[1] / dpi)
        Figure.__init__(self, figsize=figsize, dpi=dpi)
        self.patch.set_alpha(0.0)
        self.canvas = FigureCanvas(self)
        w, h = self.get_figwidth(), self.get_figheight()
        ax = self.add_axes([0, 0, 1, 1], frameon=False, autoscale_on=False, aspect='auto')
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_ylim(-1, 1)
        self.ax = ax


    def set_bbox(self, bbox):
        self.ax.set_xlim(bbox[0], bbox[2])
        self.ax.set_ylim(bbox[1], bbox[3])


    def add_patch(self, patch):
        self.ax.add_patch(patch)
        if hasattr(patch, 'text'):
            self.ax.text(patch._textx, patch._texty, patch.text, va='top', ha='left')

    def save(self, filename):
        self.canvas.print_figure(filename, dpi=self.dpi)
    savefig = save


class GeneFigure(SimpleFigure):
    def __init__(self, figsize=(512, 256), dpi=96):
        SimpleFigure.__init__(self, figsize, dpi)
        self.ax.set_ylim(-1.5, 1.5)

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

    gf.set_xlim(0, 300)
    gf.save('t.png')
