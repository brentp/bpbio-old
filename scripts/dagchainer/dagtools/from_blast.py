import sys
import gtpym
import operator

def to_order(gff):
    seqids = sorted(gff.seqids)
    order = {}
    i = 0
    for seqid in seqids:
        feats = sorted(gff.get_features_for_seqid(seqid), key=operator.attrgetter('start'))
        for feat in feats:
            order[feat.attribs["ID"]] = i
            i += 1
    return order

def blast_to_dag(blast_file, query, subject, qgff_file, sgff_file, qdups,
                 sdups, order, condense):
    """
    if order is true, we return the index of the gene, instead of the basepair position
    """
    if qdups:
        qdups = frozenset(x.strip() for x in open(qdups))
    if sdups:
        sdups = frozenset(x.strip() for x in open(sdups))

    qgff = gtpym.FeatureIndexMemory(qgff_file)
    sgff = gtpym.FeatureIndexMemory(sgff_file)

    if order:
        qorder = to_order(qgff)
        sorder = to_order(sgff)

    qorg = query   + "_"
    sorg = subject + "_"
    seen = {}
    n_qdups = 0
    n_sdups = 0
    for line in open(blast_file):
        line = line.split("\t")

        if line[0] in qdups: n_qdups += 1; continue
        if line[1] in sdups: n_sdups += 1; continue
        qfeat = qgff[line[0]]
        sfeat = sgff[line[1]]

        if condense:
            qname = qfeat.attribs["ID"]
            sname = sfeat.attribs["ID"]
            key = qname + "@" + sname
            if key in seen: continue
            seen[key] = True
        
        if order:
            qo = qorder[qfeat.attribs["ID"]]
            so = sorder[sfeat.attribs["ID"]]
            print "\t".join(map(str, [
                 qorg + qfeat.seqid, line[0], qo, qo,
                 sorg + sfeat.seqid, line[1], so, so]))

        else:
            print "\t".join(map(str, [
                 qorg + qfeat.seqid, line[0], qfeat.start, qfeat.end
                ,sorg + sfeat.seqid, line[1], sfeat.start, sfeat.end]))

    if qdups:
        print >>sys.stderr, "removed %i dups from query  " % n_qdups
        if n_qdups == 0:
            raw_input("didnt remove any dups from teh query, this is Baad!!!...  press a key")
    if sdups:
        print >>sys.stderr, "removed %i dups from subject" % n_sdups
        if n_sdups == 0:
            raw_input("didnt remove any dups from the subject, this is Baad!!!...  press a key")

def main(args):
    from optparse import OptionParser
    usage = """ 
    takes a tab-delimited blast file and converts it to the format used by
    dagchainer and tandems.py. output is to STDOUT.
    if (optional) files are given for query/subject_dups with format:
    dupa_name
    dupb_name

    .
    dupzzz_name

    then any hits containing those are removed. from the output
    """
    parser = OptionParser(usage)

    parser.add_option("-b", "--blast_file", dest="blast_file", 
                      help="the name of the blast_file", default=False)
    parser.add_option("-q", "--query",   dest="query",   
                      help="the name of the query organism")
    parser.add_option("-s", "--subject", dest="subject", 
                      help="the name of the subject organism")
    parser.add_option("-o", dest="order", action="store_true", default=False, 
                      help="if this is true, the relative positions are sent "
                      "to the dag file, rather than the basepair positions.")
    parser.add_option("--qgff",   dest="qgff", 
                      help="the path of the query gff")
    parser.add_option("--sgff", dest="sgff", 
                      help="the path of the subject subject")
    parser.add_option("--query_dups", dest="query_dups", 
                  help="file containing list of query dups", default=None)
    parser.add_option("--subject_dups", dest="subject_dups", 
                  help="file containing list of subject dups", default=None)
    parser.add_option("-c","--condense", dest="condense", 
                  help="condense duplicate blast hits", action="store_true")

    (options, _) = parser.parse_args()
    if not options.blast_file:
        sys.exit(parser.print_help())
    blast_to_dag(options.blast_file, options.query, options.subject,
                 options.qgff, options.sgff, 
                 options.query_dups, options.subject_dups, 
                 options.order, options.condense)


if __name__ == "__main__":
    if len(sys.argv) > 1:
        main(sys.argv[1:])
    else:
        main([])

