CONSTANT_MATCH_SCORE=None
MAX_MATCH_SCORE=50.0

import math
from subprocess import Popen, PIPE
import os
op = os.path
PATH = op.dirname(op.abspath(__file__))
import sys
sys.path.insert(0, PATH)
from dagtools import DagLine

def scoringF(evalue, constant_match=CONSTANT_MATCH_SCORE, max_match=MAX_MATCH_SCORE):
    if not constant_match is None:
        return constant_match

    matchScore = 10 * -math.log10(evalue);
    matchScore = int(matchScore +.5) / 10
    return max_match if matchScore > max_match else matchScore


def get_dag_line(fh):
    line = fh.readline()
    while line and line[0] == '#':
        line = fh.readline()
    if not line: return None
    return DagLine(line)

def get_name(header):
    # (reverse) Alignment #11  => 11r
    a = header.split('Alignment #')[1].split()[0]
    num_genes = int(header.split()[-1].rstrip(')'))
    return a + ('.%i' % num_genes) +  ('.r' if 'reverse' in header else '.f')


def get_meta_gene(fh, header=[None]):
    if header[0] is None:
        header[0] = fh.readline()
    line = fh.readline()
    genes = []
    while line and line[0] != "#": 
        d = DagLine(line)
        genes.append(d)
        line = fh.readline()
    if len(genes) == 0: return None
    l = header[0]
    name = get_name(header[0])
    # save the next header.
    header[0] = line

    reverse = name.endswith('r')

    a_start = min(g.a_start for g in genes)
    a_end   = max(g.a_end for g in genes)

    b_start = min(g.b_start for g in genes)
    b_end   = max(g.b_end for g in genes)
    if reverse: b_start, b_end = b_end, b_start

    d = {'a_seqid': genes[0].a_seqid,
         'b_seqid': genes[0].b_seqid,
         'a_accn': 'a' + name,
         'b_accn': 'b' + name,
         'a_start': a_start, 
         'b_start': b_start, 
         'a_end': a_end, 
         'b_end': b_end, 
         'evalue': 1e-250}
    return DagLine.from_dict(d)
    

def parse_file(dag_file, evalue_cutoff, metagene=False):

    accn_info = {}
    matches = {}
    fh = open(dag_file)
    dag = True
    while dag:
        if metagene:
            dag = get_meta_gene(fh)
            if dag is None: break

        else:
            dag = get_dag_line(fh)
            if dag is None: break

        if dag.evalue >= evalue_cutoff: continue
        if dag.a_accn == dag.b_accn: continue
        
        if not dag.a_accn in accn_info:
            mid = int((dag.a_start + dag.a_end + 0.5) / 2)
            a_feat = {'accn': dag.a_accn, 'start': dag.a_start, 'end': dag.a_end, 'mid': mid, 'seqid': dag.a_seqid}
            accn_info[dag.a_accn] = a_feat
        else:
            a_feat = accn_info[dag.a_accn]

        if not dag.b_accn in accn_info:
            mid = int((dag.b_start + dag.b_end + 0.5) / 2)
            b_feat = {'accn': dag.b_accn, 'start': dag.b_start, 'end': dag.b_end, 'mid': mid, 'seqid': dag.b_seqid}
            accn_info[dag.b_accn] = b_feat
        else:
            b_feat = accn_info[dag.b_accn]
    
        # always sort by seqid and order. 
        if dag.a_seqid > dag.b_seqid:
            a_feat, b_feat = b_feat, a_feat

        elif dag.a_seqid == dag.b_seqid and a_feat['mid'] > b_feat['mid']:
            a_feat, b_feat = b_feat, a_feat


        seqid_key = (a_feat['seqid'], b_feat['seqid'])
        if not seqid_key in matches: matches[seqid_key] = {}
        these_matches = matches[seqid_key]

        if a_feat['accn'] < b_feat['accn']:
            accn_key = a_feat['accn'], b_feat['accn']
        else:
            accn_key = b_feat['accn'], a_feat['accn']

        if accn_key in these_matches:
            if dag.evalue < these_matches[accn_key]['evalue']: these_matches[accn_key]['evalue'] = dag.evalue
        else:
            these_matches[accn_key] = {'A': a_feat, 'B': b_feat, 'evalue': dag.evalue}

    return matches

# TODO: do a filter on the final output by count of repeats and % of a dag group that
# is made up of repeats.


def run_dag_chainer(a_seqid, b_seqid, filename, matches, reverse, options,
                    child_conn,
                   dagchainer=os.path.join(os.path.abspath(os.path.dirname(__file__)), "dagchainer")):
    """
    calls dagchainer and yields groups of matches
    """
    o = options
    cmd = "%(dagchainer)s -G %(gap_length)s -O %(gap_init)s -E %(gap_extend)s -S " +\
          "%(min_score)s -D %(max_dist)s  -F %(filename)s %(reverse)s" # > %(tmp_file)s";
    cmd = cmd % dict(gap_length=o.gap_dist, gap_init=o.gap_init, 
                     gap_extend=o.gap_extend, min_score=o.min_score, 
                     max_dist=o.gap_max, filename="-", reverse=reverse,
                    dagchainer=dagchainer)
    num2pair = matches.values()
    process = Popen(cmd, stdin=PIPE, stdout=PIPE, bufsize=8*4096, shell=True)
    write = process.stdin.write
    for i, pair in enumerate(num2pair):
        write("%i\t%i\t%i\t%.4f\n" % (i, pair['A']['mid'], pair['B']['mid'], scoringF(pair['evalue'])))
    process.stdin.close()

    header = None
    all_data = [] # added instead of yield to allow parallelization.
    data = []
    for line in process.stdout:
        if line[0] == ">":
            if header is None:
                header = line[1:].strip()
            else:
                if len(data) >= o.min_aligned_pairs:
                    #yield header, data
                    all_data.append((header, data))
                header = line[1:].strip()
                data = []
            continue

        index, pair_id, pos1, pos2, match_score, dag_chain_score = line.strip().split()
        pair = num2pair[int(pair_id)]
        data.append({'pair': pair, 'dag_score': float(dag_chain_score)})

    if len(data) >= o.min_aligned_pairs:
        #yield header, data
        all_data.append((header, data))
    child_conn.send(all_data)
    child_conn.close()
    
def print_alignment(header, group, opts):
    header_fmt = "## alignment %s vs. %s %s (num aligned pairs: %i)"

    d = group[0]['pair']
    print header_fmt % (d['A']['seqid'], d['B']['seqid'], header, len(group))

    for pair_dict in group:
        A = pair_dict['pair']['A']
        B = pair_dict['pair']['B']
        print "%s\t%s\t%d\t%d\t%s\t%s\t%d\t%d\t%e\t%d" % (\
                     A['seqid'], A['accn'], A['start'], A['end'],
                     B['seqid'], B['accn'], B['start'], B['end'],
                     pair_dict['pair']['evalue'], pair_dict['dag_score'])

if __name__ == "__main__":
    import optparse, sys
    p = optparse.OptionParser()

    p.add_option('-i', dest='dag', help="""dag file with format
a_seqid<tab>a_accn<tab>a_start<tab>a_end<tab>b_seqid<tab>b_accn<tab>b_start<tab>b_end<tab>e-value""")
    p.add_option('-o', dest='gap_init', type='float', default=0, 
                help="gap open penalty")
    p.add_option('-e', dest='gap_extend', type='float', default=-3, 
                help="gap extension penalty")

    p.add_option('-x', dest="min_score", type='float', default=None,
                 help="minimum alignment score. alignment stops when score " + \
                 " below this value")
    p.add_option('-g', dest='gap_dist', type='float', default=10000, 
                help="averge distance expected between 2 syntenic genes")

    p.add_option('-D', dest='gap_max', type='float', default=200000,
                help="maximum distance between 2 matches in basepairs")
    p.add_option('-E', dest='evalue', type='float', default=1e-5,
                help="maximum evalue.")
    p.add_option('-A', dest='min_aligned_pairs', type='int', default=6,
                help="minimum number of pairs to be considered a diagonal")

    p.add_option('-I', dest='ignore_tandem', default=False, action='store_true',
                help="ignore tandems in self-self comparisons.")

    p.add_option('-M', dest='max_match_score', type='float', default=50,
                help="maximum score to be assigned to a match")
    p.add_option('--meta_genes', dest='meta_genes', default=False, 
                 action='store_true', help="""\
                 if this flag is set, then the input is assumed to be
                 a file that is normally output by dagchainer, containing
                 genes grouped by diagonal. each group is merged into a single
                 'meta' gene and the resulting 'meta'-genes are run as normal
                 through dagchainer.""")

    opts, _ = p.parse_args() 

    if not opts.dag:
        sys.exit(p.print_help())

    if opts.min_score is None:
        opts.min_score = int(opts.min_aligned_pairs * 0.5 * opts.max_match_score)

    all_matches = parse_file(opts.dag, opts.evalue, opts.meta_genes)

    from multiprocessing import Process, Pipe as mPipe
    
    filename = "-"
    for (a_seqid, b_seqid), matches in sorted(all_matches.iteritems()):

        parent_connf, child_connf = mPipe()
        pf = Process(target=run_dag_chainer, args=(a_seqid, b_seqid, filename, matches, "", opts, child_connf))
        pf.start()

        parent_connr, child_connr = mPipe()
        pr = Process(target=run_dag_chainer, args=(a_seqid, b_seqid, filename, matches, "-r", opts, child_connr))
        pr.start()

        for header, group in parent_connf.recv():
            print_alignment(header, group, opts)
        for header, group in parent_connr.recv():
            print_alignment("(reverse) " + header, group, opts)

        pr.join()
        pf.join()
