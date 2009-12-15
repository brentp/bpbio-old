#!/usr/bin/python
import sys
sys.stdout = sys.stderr
from cgi import parse_qsl
import tables
import os
os.environ[ 'HOME' ] = '/tmp/'
import matplotlib
matplotlib.interactive(0)
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from cStringIO import StringIO

import glob
DIR = os.path.abspath(os.path.dirname(__file__))
f = glob.glob(DIR + "/*.h5")[0]
h5 = tables.openFile(f)

def application(env, start_response):
    start_response("200 OK", [("Content-type", "image/png")])

    p = dict(parse_qsl(env['QUERY_STRING']))


    xmin = max(0, int(p['xmin']))
    xmax = int(p['xmax'])

    chrn = "c%i" % int(p['chr'])
    vals = h5.getNode('/count/%s' % (chrn,))[xmin:xmax+1]
    ymax = 200

    io = StringIO()


    f = Figure(frameon=False)
    f.canvas = FigureCanvas(f)
    dpi = 128.
    f.set_size_inches(int(p['width'])/dpi, 300/dpi)
    ax = f.add_axes((0,0,1,1), frameon=False, xticks=(), yticks=())
    ax.set_autoscale_on(0)
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(0, ymax)
    ax.plot(range(xmin, xmax+1), vals)
    ax.plot([xmin,xmax],[50,50],'r')
    f.savefig(io, format='png', dpi=dpi)
    io.seek(0)

    return [io.read()]
