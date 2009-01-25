import os
from subprocess import Popen
import sys

import logging
logging.basicConfig(level=logging.DEBUG)
log = logging.getLogger(os.path.basename(__file__))

def get_blast_file(qfasta, sfasta, out_dir=None):
    q = os.path.basename(qfasta)
    s = os.path.basename(sfasta)
    blast_file = q[:q.rfind(".")] + "_vs_" + s[:s.rfind(".")] + ".blast"
    if not out_dir: return blast_file
    if out_dir is True:
        d = os.path.dirname(qfasta)
        return os.path.join(d, blast_file)
    return os.path.join(out_dir, blast_file)

def is_current_file(a, b):
    if not (os.path.exists(a) and os.path.exists(b)): return False
    am = os.stat(a).st_mtime
    bm = os.stat(b).st_mtime
    return am > bm

def is_same_blast_params(blast_file, cmd):
    """ when using the blast() below, a .cmd file is
    written, this function checks that file to see if the
    current command is the same as that. if so, and the fasta
    files are up to date, the blast is not done as it's up to date"""
    params_file = blast_file + ".cmd"
    if not os.path.exists(params_file): return False
    return cmd.strip() == open(params_file).read().strip()

def sh(cmd):
    """ run a commmand in the shell"""
    log.debug(cmd)
    proc = Popen(cmd, stdout=sys.stdout, stderr=sys.stderr, shell=True)
    r = proc.communicate()

def add_dash(params):
    d = {}
    for k, v in params.items():
        if k.startswith("-"):
            d[k] = v
        elif len(k) > 1:
            d["--" + k] = v
        else:
            d["-" + k] = v
    return d


def is_protein_db(blast_cfg):
    return blast_cfg["p"] in ("blastx", "blastp")
   

def blast(_blast_cfg, blastall="/usr/bin/blastall", full_name=False):
    blast_cfg = _blast_cfg.copy()
    check_args(blast_cfg)
    q_fasta = blast_cfg["i"]
    s_fasta = blast_cfg["d"]

    format_db = blastall.replace("blastall", "formatdb")

    protein = "T" if is_protein_db(blast_cfg) else "F"
    cmd = "%(format_db)s -i %(s_fasta)s -p %(protein)s" % locals()

    if not is_current_file(s_fasta + ".nin", s_fasta):
        sh(cmd)
    else:
        log.warn("NOT running cmd:\n%s\n because %s.nin is up to date" % (cmd, s_fasta))
    blast_file = ""
    to_query_dir = blast_cfg.get("o", "F").upper() != "F"
    if blast_cfg.get("o", "F").upper() not in ("T", "F"):
        to_query_dir = blast_cfg["o"]
    blast_file = get_blast_file(q_fasta, s_fasta, to_query_dir)

    if full_name:
        blast_file = blast_file.rstrip(".blast") \
              + "_params__" \
              + "__".join(["%s_%s" % p for p in sorted(blast_cfg.items())
                           if not p[0] in ("i", "d")]) \
              + ".blast"
    blast_cfg.update({"o": blast_file})

    params = add_dash(blast_cfg)

    params = ["%s %s" % (p, v) for p, v in sorted(params.items())]
    cmd = blastall + " " + " ".join(params)

    if not (is_current_file(blast_file, q_fasta) \
                and is_current_file(blast_file, s_fasta) \
                and is_same_blast_params(blast_file, cmd)):

        fh = open(blast_file + ".cmd", "w")
        fh.write(cmd)
        fh.close()
        sh(cmd)
        if os.path.exists(blast_file):
            lines = sum(1 for line in open(blast_file))
            log.debug("\n\n%s lines of blast output sent to %s" % (lines, blast_file))
        else:
            log.errorr("\n\nERROR: blast not run")
    else:
        log.error("NOT running cmd:\n%s\n because %s is up to date" % (cmd, blast_file))

def check_args(args):
    if not "p" in args: args["p"] = "blastn"
    assert "i" in args, "need to specify a query fasta"
    assert "d" in args, "need to specify a query fasta"
    if not "a" in args: args["a"] = "4"


if __name__ == "__main__":
    if len(sys.argv) == 1 or sys.argv[1] in ("-h", "--help"):
        sh("/usr/bin/blastall")
        print """\
   this script will generally do exactly the same as blastall
   except it will attempt to create the name of the output blast
   file from the input fasta files, and the blast parameters.
   it will also run formatdb with the correct -p parameter based
   on the type of blast requested.
   it also saves a file: a_vs_b.blast.cmd that stores the exact
   blast command used to generate a_vs_b.blast
   in addition, if a.fasta an b.fasta are older than a_vs_b.blast
   and the params have not changed, the blast will not run as
   the current blast file is up to date.

   additional args provided by this script are:
            -o T
         or 
            -o F

         in which case the blast output file will be created
         from the names of the input fasta files and placed in the
         directory of the query fasta
         in the latter case, the blast file will go to the current
         direcotry

            --full_name T

         if specified, this will include the blast params in the name
         of the output blast file. e.g.: a__vs_b__params__m_8__W_15.blast
           for -m -8 -W 15
    """
        sys.exit()
    args = dict((sys.argv[i].lstrip("-") , sys.argv[i + 1].rstrip()) \
                          for i in range(1, len(sys.argv), 2))

    try:
        f = args.pop("full_name")
        full_name = not f.lower() in ("f", "0")
    except:
        full_name = False

    if not "i" in args:
        print "need to specify a query fasta (-i)"
        sys.exit()
    if not "d" in args:
        print "need to specify a subject fasta (-d)"
        sys.exit()
    
    blast(args, full_name=full_name)
