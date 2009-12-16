CONSTANT_MATCH_SCORE=None
MAX_MATCH_SCORE=50.0

import math
from subprocess import Popen, PIPE
import os


def scoringF(evalue, constant_match=CONSTANT_MATCH_SCORE, max_match=MAX_MATCH_SCORE):
    if not constant_match is None:
        return constant_match

    matchScore = 10 * -math.log10(evalue);
    matchScore = int(matchScore +.5) / 10
    return max_match if matchScore > max_match else matchScore

def parse_file(dag_file, evalue_cutoff):

    dag_cols = ('a_seqid', 'a_accn', 'a_start', 'a_end', 'b_seqid', 'b_accn', 'b_start', 'b_end', 'evalue')
    accn_info = {}
    matches = {}
    for line in open(dag_file):
        if line[0] == '#': continue
        tline = line[:-1].split("\t")

        tline[2:4] = map(int, tline[2:4])
        tline[6:8] = map(int, tline[6:8])
        tline[8] = float(tline[8])
        dag = dict(zip(dag_cols, tline))
        if dag['evalue'] >= evalue_cutoff: continue
        if dag['a_accn'] == dag['b_accn']: continue
        if dag['evalue'] < 1e-250: dag['evalue'] = 1e-250
        
        if not dag['a_accn'] in accn_info:
            mid = int((dag['a_start'] + dag['a_end'] + 0.5) / 2)
            a_feat = {'accn': dag['a_accn'], 'start': dag['a_start'], 'end': dag['a_end'], 'mid': mid, 'seqid': dag['a_seqid']}
            accn_info[dag['a_accn']] = a_feat
        else:
            a_feat = accn_info[dag['a_accn']]

        if not dag['b_accn'] in accn_info:
            mid = int((dag['b_start'] + dag['b_end'] + 0.5) / 2)
            b_feat = {'accn': dag['b_accn'], 'start': dag['b_start'], 'end': dag['b_end'], 'mid': mid, 'seqid': dag['b_seqid']}
            accn_info[dag['b_accn']] = b_feat
        else:
            b_feat = accn_info[dag['b_accn']]
    
        # always sort by seqid and order. 
        if dag['a_seqid'] > dag['b_seqid']:
            a_feat, b_feat = b_feat, a_feat

        elif dag['a_seqid'] == dag['b_seqid'] and a_feat['mid'] > b_feat['mid']:
            a_feat, b_feat = b_feat, a_feat


        seqid_key = (a_feat['seqid'], b_feat['seqid'])
        if not seqid_key in matches: matches[seqid_key] = {}
        these_matches = matches[seqid_key]

        accn_key = tuple(sorted([a_feat['accn'], b_feat['accn']]))

        if accn_key in these_matches:
            if dag['evalue'] < these_matches[accn_key]['evalue']: these_matches[accn_key]['evalue'] = dag['evalue']
        else:
            these_matches[accn_key] = {'A': a_feat, 'B': b_feat, 'evalue': dag['evalue']}

    return matches

TMPDIR='/tmp/'
def gen_matches_by_seqid(matches):
    for a_seqid, b_seqid in sorted(matches):
        these_matches = matches[(a_seqid, b_seqid)]
        # TODO: redundant return values and dict. but need to keep them ordered.
        # TODO: make filename a param if needed. (or just wirte to file if debug is on)
        yield a_seqid, b_seqid, "-", these_matches


def run_dag_chainer(a_seqid, b_seqid, filename, matches, reverse, options,
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
    #print >>sys.stderr, cmd
    num2pair = matches.values()
    
    process = Popen(cmd, stdin=PIPE, stdout=PIPE, shell=True)
    for i, pair in enumerate(num2pair):
        print >>process.stdin, "%i\t%i\t%i\t%.4f"% (i, pair['A']['mid'], pair['B']['mid'], scoringF(pair['evalue']))
    process.stdin.close()

    header = None
    data = []
    for line in process.stdout:
        if line[0] == ">":
            if header is None:
                header = line[1:].strip()
            else:
                if len(data) >= o.min_aligned_pairs:
                    yield header, data
                header = line[1:].strip()
                data = []
            continue

        index, pair_id, pos1, pos2, match_score, dag_chain_score = line.strip().split()
        pair = num2pair[int(pair_id)]
        data.append({'pair': pair, 'dag_score': float(dag_chain_score)})

    if len(data) >= o.min_aligned_pairs:
        yield header, data

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

    opts, _ = p.parse_args() 

    if not opts.dag:
        sys.exit(p.print_help())

    if opts.min_score is None:
        opts.min_score = int(opts.min_aligned_pairs * 0.5 * opts.max_match_score)

    all_matches = parse_file(opts.dag, opts.evalue)
    for match_info in gen_matches_by_seqid(all_matches):
        a_seqid, b_seqid, filename, matches = match_info
        #TODO: forward and reverse.
        for header, group in run_dag_chainer(a_seqid, b_seqid, filename, matches, "", opts):
            print_alignment(header, group, opts)
        for header, group in run_dag_chainer(a_seqid, b_seqid, filename, matches, "-r", opts):
            print_alignment("(reverse) " + header, group, opts)
