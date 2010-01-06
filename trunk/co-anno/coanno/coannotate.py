import simplejson
import logging
import bblast
import collections
from gtpym import FeatureIndexMemory
import operator
import os
import sys
import string
from biostuff import BlastLine
from pyfasta import Fasta
import mask_features as mf


def gff_feat_to_str(f):
    r"""
    >>> f = dict(name='name', start=2, stop=3, strand='+', chr='chr3', 
    ...           type='CDS', attrs={'match':'amatch', 'ID': 'at2g26540'})
    >>> gff_feat_to_str(f)
    'chr3\tucb\tCDS\t2\t3\t.\t+\t.\tID=at2g26540;match=amatch'
    """
    attrs = ";".join(["%s=%s" % (k,v) for k, v in f['attrs'].iteritems()])
    strand = f['strand'] in (1, '1', '+') and '+' or '-'
    return "\t".join(map(str, (f['chr'], 'ucb', f.get('type', 'gene'), f['start'], f['stop'], ".", strand, ".", attrs)))


logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger('coanno')


MATCH_JOIN = "#"

def run_mask(cfg):
    for g in ("genome_a", "genome_b"):
        mf.main(cfg[g]["gff"], cfg[g]["fasta"], cfg["default"]["mask_order"])
        log.debug("MASKED: %s" % cfg[g]['name'])

def update_cfg(cfg):
    """
    set some defaults.

    >>> cfg = {'default': {}}
    >>> update_cfg(cfg)
    >>> cfg['default']['min_overlap']
    30
    """

    if not "overlap" in cfg["default"]:
        cfg["default"]["min_overlap"] = MIN_OVERLAP
    if not "coverage" in cfg["default"]:
        cfg["default"]["min_pct_coverage"] = MIN_PCT_COVERAGE


def main(cfg_path):
    cfg = simplejson.load(open(cfg_path))
    update_cfg(cfg)

    run_mask(cfg)

    run_blasts(cfg)

    dispatch(cfg)
    if cfg["default"].get("reciprocal"):
        dispatch(cfg, flip=True)

MIN_OVERLAP = 30
MIN_PCT_COVERAGE = 0.40

def merge_overlapping(new_genes, min_overlap, min_pct_coverage):
    """
    require at least `overlap` basepairs worth of overlap
    >>> list(merge_overlapping({}, 30, 0.4))
    []

    """

    merged = 0
    # since we have to back track sometimes
    # keep track of stuff on other strand.
    #dbg = open("genesets.txt", "w")

    for achr, genes in sorted(new_genes.iteritems()):
        genes = sorted([g for _name, g in genes.items()], key=operator.itemgetter('start'))
        i, j = 0, 1
        used_genes = set([])
        while i < len(genes):

            # it was probably seen on the other strand.
            if i in used_genes: i += 1; j = i + 1; continue
            g = genes[i].copy()
            matches = g['attrs']['match']
            new_stop = g['stop']
            #print >>dbg, g

            last_gene = None # keep track of gene on other strand.

            #  PARAMETER: require 30bp overlap.
            while j < len(genes) and g['stop'] > genes[j]['start'] + min_overlap:
                # can only merge genes with genes, and cds with cds.
                if j in used_genes or genes[j]['type'] != g['type']: j += 1; continue
                #print g['name'], genes[j]['name']
                minstart = g['start']
                maxstop = max(new_stop, genes[j]['stop'])
                maxstart = genes[j]['start']
                minstop = min(new_stop, genes[j]['stop'])
                coverage = (minstop - maxstart) / (1.0 * maxstop - minstart)
                if genes[j]['strand'] != g['strand'] and coverage < min_pct_coverage:
                    last_gene = j
                    j += 1
                    continue
                new_stop = maxstop
                #print >>dbg, genes[j]
                used_genes.update([j])
                # blast doesnt like looong names.
                if len(matches) < 64:
                    matches += MATCH_JOIN + genes[j]['attrs']['match']
                j += 1
                merged += 1
            #print >>dbg,""
          
            g['stop'] = new_stop
            g['attrs']['match'] = matches[:64]
            fix_name(g)
            assert len(g['name']) <= 80, (len(g['name']), g)


            yield g
            used_genes.update([i])
            i = last_gene or j
            j = i + 1
    print >>sys.stderr, "merged %i genes" % merged

def fix_name(g):
    """the name is like: org_chr_start_stop but we may have changed start, stop

    >>> f = dict(name='name', start=2, stop=3, strand='+', chr='chr3', 
    ...           type='CDS', attrs={'match':'amatch', 'ID': 'at2g26540_4_5'})
    >>> fix_name(f)
    >>> f
    {'start': 2, 'chr': 'chr3', 'name': 'at2g26540_2_3', 'type': 'CDS', 'stop': 3, 'strand': '+', 'attrs': {'rname': 'at2g26540_2_3', 'ID': 'at2g26540_2_3', 'match': 'amatch'}}
    
    `"""
    new_name = g['attrs']['ID']
    new_name = new_name.rstrip(string.digits).rstrip("_").rstrip(string.digits) + "%i_%i"
    new_name %= (g['start'], g['stop'])
    g['name'] = g['attrs']['rname'] = g['attrs']['ID'] = new_name

def dispatch(cfg, flip=False):
    """
    >>> cfg = simplejson.loads(open("tests/data/athaliana_athaliana.json").read())
    >>> update_cfg(cfg)
    >>> dispatch(cfg)

    """
    if flip is False:
        akey, bkey = "genome_a", "genome_b"
    else:
        akey, bkey = "genome_b", "genome_a"

    afasta = cfg[akey]["fasta"]
    bfasta = cfg[bkey]["fasta"]
    odir = cfg["default"]["out_dir"]
    assert os.path.exists(odir), "Need to create directory: %s" % odir
    min_len = cfg["default"]["min_len"]
    min_pct_coverage = cfg["default"]["min_pct_coverage"]

    out_gff = "missed_%s_from_%s.gff" % (cfg[bkey]["name"], cfg[akey]["name"])
    out_gff = os.path.join(cfg["default"]["out_dir"], out_gff)
    out_fh = open(out_gff, "w")
    

    a, b = fastas_for_features_vs_masked_genomic(afasta, bfasta)
    a_bnon_blast = bblast.get_blast_file(a, b, odir)
    a_b_blast = bblast.get_blast_file(afasta.replace(".fa", ".features.fa"),
                                      bfasta.replace(".fa", ".features.fa"), odir)
    i = 0
    new_genes = collections.defaultdict(dict)
    for new_gene in find_missed(cfg[bkey]["name"],
                            cfg[akey]["gff"], cfg[bkey]["gff"],
                            a_bnon_blast, a_b_blast, min_pct_coverage, flip=flip, min_len=min_len):
        #if new_gene['name'] in seen: continue
        n, achr = new_gene['name'], new_gene['chr']
        # found it 2x, add the match to the list.
        if n in new_genes[achr]:
            new_genes[achr][n]['attrs']['match'] += MATCH_JOIN + new_gene['attrs']['match']
        else:
            new_genes[achr][n] = new_gene
    # TODO: merge overlapping. only need to check chr, start, stops. so
    # sortable without index.
    merged_genes = merge_overlapping(new_genes, cfg['default']['min_overlap'],
                                     cfg['default']['min_pct_coverage'])

    merged_genes = exclude_genes_in_high_repeat_areas(merged_genes, bfasta)
    print >>out_fh, "##gff-version 3"
    for new_gene in merged_genes:
        print >>out_fh, gff_feat_to_str(new_gene)
        i += 1
    log.debug("created %i new features in %s" % (i, out_gff))


def exclude_genes_in_high_repeat_areas(merged_genes, bfasta):
    #print "FASTA:", afasta

    f = Fasta(bfasta)
    skipped = 0
    for gene in merged_genes:
        # get the total sequence length.
        seq = str(f[gene['chr']][gene['start']:gene['stop']])
        tot = len(seq)
        # and the lenght of real sequence.
        seq = seq.upper().replace('N', '').replace('X','')
        # if it's not > 80% good sequence, just skip it.
        if float(len(seq))/tot < .85: skipped+=1; continue
        yield gene
    log.info("removed %i (otherwise) new genes in masked areas" % skipped)


def dist(a, b):
    """
    >>> class O():
    ...    def __init__(self, start, stop):
    ...        self.sstart, self.sstop = start, stop

    >>> a, b, c = O(12, 15), O(10, 11), O(16, 17)
    >>> dist(a, b)
    1
    >>> dist(a, c)
    1
    >>> dist(b, c)
    5
    """

    d1 = abs(a.sstart - b.sstart)
    d2 = abs(a.sstop - b.sstart)
    d3 = abs(a.sstop - b.sstop)
    d4 = abs(a.sstart - b.sstop)
    return min(d1, d2, d3, d4)

def partition(slist, max_dist=1000):
    """
    return a list of lists where the hits are sorted by location and things 
    farther than max_dist apart are in separate lists.

    >>> class O():
    ...    def __init__(self, start, stop):
    ...        self.sstart, self.sstop = start, stop
    ...    def __repr__(self): return "O(%i, %i)" % (self.sstart, self.sstop)

    >>> partition([1])
    [[1]]

    >>> slist = [O(2, 5), O(1, 10), O(6, 7), O(160, 170), O(171, 180)]
    >>> partition(slist, max_dist=100)
    [[O(1, 10), O(2, 5), O(6, 7)], [O(160, 170), O(171, 180)]]

    >>> slist = [O(2, 5), O(1, 10), O(6, 7), O(26, 27), O(160, 170), O(171, 180)]
    >>> partition(slist, max_dist=10)
    [[O(1, 10), O(2, 5), O(6, 7)], [O(26, 27)], [O(160, 170), O(171, 180)]]

    """
    if not isinstance(slist, list):
        slist = list(slist)
    if len(slist) == 1: return [slist]
    slist.sort(key=operator.attrgetter('sstart'))
    lists = []
    current_list = [slist[0]]
    for a in slist[1:]:
        for b in current_list:
            if dist(a, b) < max_dist:
                current_list.append(a)
                break
        else:
            lists.append(current_list[:])
            current_list = [a]
    lists.append(current_list)
    return lists

def find_missed(sorg, qgff_file, sgff_file, q_snon_blast, q_s_blast,
                min_pct_coverage, min_len=30, flip=False):
    """ e.g.:
        >>> find_missed("papaya", "grape.gff", "papaya.gff",
        ...         "grape.features_vs_papaya.genomic.masked.blast",
        ...         "grape.features_vs_papaya.features.blast", 0.30) # doctest: +ELLIPSIS
        <generator object find_missed at ...>

    to find papaya mised exons. and:
        >>> find_missed("grape", "papaya.gff", "grape.gff",
        ...         "papaya.features_vs_grape.genomic.masked.blast",
        ...         "grape.features_vs_papaya.features.blast", 0.30, flip=True) # doctest: +ELLIPSIS
        <generator object find_missed at ...>

    to find grape missed exons.

    min coverage means if it finds a bunch of little spread out hits, it
    discards them unless they cover at least (for example) 0.4 == 40% of the
    total area. so hits like [105, 120], [190, 205] only cover 
       15 + 15 / (205 - 105) = 30%
    """
    sgff = FeatureIndexMemory(sgff_file)
    qgff = FeatureIndexMemory(qgff_file)

    # grouped_by_q has all the subject genomic hits mapped to the query feature.
    # and the subject start, stop are chromosomal positions.
    name = sorg + "_%(chr)s_%(start)i_%(stop)i"
    grouped_by_q = grouper(q_snon_blast)
    #print len(grouped_by_q)


    for qname, sdict in sorted(grouped_by_q.iteritems()):
        try:
            qfeat = qgff[qname]
        except:
            print >>sys.stderr, "\n", qname, sdict, qgff_file, sgff_file, "\n"
            raise
        qstrand = qfeat.strand == '-' and -1 or 1


        for schr, big_slist in sorted(sdict.iteritems()):
            slists = partition(big_slist)

            cover = 0.0
            for slist in slists:
                sstrand = slist[0].sstart < slist[0].sstop and 1 or -1
                if sstrand == 1:
                    cover += sum(x.sstop - x.sstart for x in slist)
                    sstart = min([x.sstart for x in slist])
                    sstop  = max([x.sstop for x in slist])
                else:
                    sstart = min([x.sstop for x in slist])
                    sstop  = max([x.sstart for x in slist])
                    cover += sum(x.sstart - x.sstop for x in slist)
                if abs(sstop - sstart) < min_len: continue

                # the hsps have to not be toooooo sparse
                if cover / (sstop - sstart) < min_pct_coverage: continue
                sstrand *= qstrand
                sname = name % dict(chr=schr, start=sstart, stop=sstop)



                feat = dict(name=sname , start=sstart , stop=sstop
                           , chr=schr , type="gene"
                           , strand= sstrand== 1 and "+" or "-"
                           , attrs={"match": qfeat.get_attribute('ID'), "ID": sname,
                               'rname':sname})

                # check if it's inside an existing gene. in which case, 
                # call it a cds and give the the same name as the parent.
                parent = [s for s in sgff.get_features_for_range(feat['start'],
                    feat['stop'], feat['chr'])] # if s.strand == feat['strand']]

                if len(parent) == 0:
                    yield feat
                else:
                    # it overlaps 2 of them, in between 2 genes. what to do?
                    # just use the first 1 for now.
                    #assert len(parent) == 1, (parent, parent[0].attribs["ID"],
                    #        parent[1].attribs["ID"])
                    parent = parent[0]
                    # when doing self-self dont want to add an annotation based
                    # on the same gene.
                    if parent.get_attribute("ID") == qfeat.get_attribute("ID"):
                        continue
                    feat['name'] = feat['attrs']['rname'] = feat["attrs"]["Parent"] = parent.get_attribute("ID")
                    feat['type'] = 'CDS'
                    feat['strand'] = parent.strand
                    #del feat['attrs']['match']
                    feat['attrs']['anno'] = "ME" 
                    yield feat

                

def grouper(blast_file):
    """\
    group all subjects to a single query. so for
        grape.features_vs_papaya.genomic.masked.blast
    group all the papaya hits to the grape query"""
    g = collections.defaultdict(dict)
    for sline in open(blast_file):
        b = BlastLine(sline)
        # this removes low-copy, transposons (length > 200, percent_id > 98)
        if b.pctid > 98.0 and b.hitlen > 200: continue
        if not b.subject in g[b.query]: g[b.query][b.subject] = []
        g[b.query][b.subject].append(b)
    return g

def fastas_for_features_vs_masked_genomic(afasta, bfasta):
    """
    >>> fastas_for_features_vs_masked_genomic("a.fasta", "b.fasta")
    ('a.features.fasta', 'b.genomic.masked.fasta')
    """
    a = afasta.replace(".fa", ".features.fa")
    b = bfasta.replace(".fa", ".genomic.masked.fa")
    return a, b

def run_blasts(config, test=False):
    """
    >>> cfg = simplejson.loads(open("tests/data/athaliana_athaliana.json").read())
    >>> run_blasts(cfg, test=True)

    """
    blast_cfg = config["blast"]
    blast_cfg.update({"m": 8, "p": "blastn"})

    out_dir = config["default"]["out_dir"]
    blast_cfg["o"] = out_dir
    blastall = config["default"]["blast_path"]
    # first blast a features to b genomic with features masked.
    agenes, bfasta =\
                fastas_for_features_vs_masked_genomic(config["genome_a"]["fasta"],
                                                      config["genome_b"]["fasta"])
    blast_cfg["i"] = agenes
    blast_cfg["d"] = bfasta

    blast_log = config["default"].get("blast_log")
    if blast_log:
        blast_log = config["genome_a"]["name"] + "_" \
                + config["genome_b"]["name"] + "_blast.log"

    if not test:
        bblast.blast(blast_cfg, blastall=blastall, blast_log=blast_log)

    bgenes, afasta =\
                fastas_for_features_vs_masked_genomic(config["genome_b"]["fasta"],
                                                      config["genome_a"]["fasta"])

    # then vice-versa
    if config["default"].get("reciprocal"):
        blast_cfg["i"] = bgenes
        blast_cfg["d"] = afasta
        if blast_log:
            blast_log = blast_log.replace('.log', '.reciprocal.log')

        if not test:
            bblast.blast(blast_cfg, blastall=blastall, blast_log=blast_log)
    

    # then genes to genes
    #blast_cfg["i"] = agenes
    #blast_cfg["d"] = bgenes
    #bblast.blast(blast_cfg, blastall=blastall)
    


if __name__ == "__main__":
    main(sys.argv[1])
