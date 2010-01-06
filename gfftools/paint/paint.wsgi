#!/usr/bin/python

import web
import os, sys
from hashlib import md5
from web.contrib.template import render_mako
os.environ['HOME'] = '/tmp/'
import matplotlib
matplotlib.use('Agg')
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import matplotlib.pyplot as plt
from cStringIO import StringIO
#sys.stdout = sys.stderr

import to_gff
import paint2 as paint
import gt

urls = (
    '/png/(.+)', 'png',
    '/(.+)?', 'index', 
)

op = os.path
render = render_mako(directories=[op.join(op.dirname(__file__),'templates')])

SAVE_PATH = '/tmp/paint/'


def save_file(contents, hashed):
    f = os.path.join(SAVE_PATH, hashed + ".gff")
    gff_line = None
    if not os.path.exists(f):
        fh = open(f, 'wb')

        fh.write(to_gff.GFF_VERSION_HEADER + '\n')
        for gff_line in to_gff.to_gff_lines(StringIO(contents)):
            #print gff_line
            fh.write(gff_line + '\n')
        fh.close()
    else:
        gff_line = to_gff.to_gff_lines(StringIO(contents)).next()
    fchr = gff_line.split("\t")[0]
    return fchr

def get_png(hash, seqid, w, h):
    return "%s_%s_%.2f_%.2f.png" % (hash, seqid, w, h)

class png(object):
    def GET(self, pngp):
        web.header('Content-type', 'image/png')
        #w = float(web.input(width=2.2).width)
        #h = float(web.input(height=8).height)
        #seqid = web.input().seqid
        #pngpath = get_png(hashed, seqid, w, h)
        return open(op.join(SAVE_PATH, pngp), 'r').read()


class index(object):
    def GET(self, hashed=None):
        web.header('Content-type', 'text/html')

        if not hashed: return render.upload()
        if hashed.startswith('fav'): return []
        w = float(web.input(width=2.2).width)
        h = float(web.input(height=8).height)

        gff = gt.FeatureIndexMemory()
        path = op.join(SAVE_PATH, hashed + ".gff")
        web.debug("about to load gff3file:%s" % path)
        assert op.exists(path)
        gff.add_gff3file(path)
        web.debug("loaded gff3file")
        seqids = gff.get_seqids()
        seqid = web.input(seqid=seqids[0]).seqid
        pngpath = get_png(hashed, seqid, w, h)

        if not op.exists(op.join(SAVE_PATH, pngpath)):
            feats = gff.get_features_for_seqid(str(seqid))
            fig = plt.figure(figsize=(w, h))
            paint.paint_features(feats, fig)
            web.debug('done_painting')
            plt.savefig(op.join(SAVE_PATH, pngpath))

        return render.index(aseqid=seqid, hashed=hashed, seqids=seqids, pngpath="/png/" + pngpath, home=web.ctx.realhome + '/',
                width=w, height=h)
    
    def POST(self, unused=None):
        x = web.input(afile={})
        contents = x['afile'].value.strip()
        fname = x['afile'].filename

        hashed = md5(contents).hexdigest()
        fchr = save_file(contents, hashed)
        return web.seeother('/%s?seqid=%s' % (hashed, fchr))


app = web.application(urls, globals(), autoreload=True)

application = app.wsgifunc()

if __name__ == "__main__":
    app.run()



