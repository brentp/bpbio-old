"""
>>> data = \"\"\"Name\\tTAIR\\tChr\\tStart\\tStop\\tStrand
... miR156a\\tAT2G25095\\t2\\t10683551\\t10683632\\t-1
... miR156b\\tAT4G30972\\t4\\t15074951\\t15075030\\t1
... miR156c\\t\\t4\\t15413325\\t15415879\\t-1
... miR156d\\tAT5G10946\\t5\\t3456211\\t3457245\\t-1
... miR156e\\tAT5G11977\\t5\\t3867214\\t3867309\\t1\"\"\"


>>> from cStringIO import StringIO
>>> txt_data = StringIO(data)

>>> gffiter = to_gff_lines(txt_data, as_dict=True)
>>> gffiter.next()
{'start': 10683551, 'seqid': '2', 'end': 10683632, 'attrs': {'tair': 'AT2G25095', 'ID': 'miR156a'}, 'strand': '-'}

>>> arr = data.split("\\n")
>>> header_line = arr[0]
>>> sep = get_sep(header_line)
>>> sep
'\\t'

>>> h = get_header_map(header_line, sep)
>>> h
{0: 'id', 1: 'tair', 2: 'seqid', 3: 'start', 4: 'end', 5: 'strand'}

>>> l = line_to_dict(arr[1], sep, h)
>>> l
{'start': 10683551, 'seqid': '2', 'end': 10683632, 'attrs': {'tair': 'AT2G25095', 'ID': 'miR156a'}, 'strand': '-'}

>>> feature_dict_to_gff(l).split("\\t")
['2', '.', '.', '10683551', '10683632', '.', '-', '.', 'ID=miR156a;tair=AT2G25095']

"""


import sys
def get_sep(line):
    c = line.count(",")
    t = line.count("\t")
    s = line.count(" ")
    m = max(c, t, s)
    for char, sep in zip((",", "\t", " ") ,(c, t, s)):
        if m == sep: return char
    raise Exception("cant determine separator")

def get_header_map(line, sep):
    line = line.rstrip().lower().split(sep)
    assert 'start' in line
    assert 'stop' in line or 'end' in line
    assert 'seqid' in line or 'chr' in line or 'chromosome' in line
    line = (l if l != 'stop' else 'end' for l in line)
    line = (l if (l != 'chr' and l != 'chromosome') else 'seqid' for l in line)
    line = (l if l != 'name' else 'id' for l in line)
    return dict(enumerate(line))

GFF_VERSION_HEADER = '##gff-version 3'
GFF_HEADERS = {'seqid': str, 'source': str, 'type': str, 'start': int, 'end':
               int, 'score':float, 'strand': str, 'phase': str}
GFF_ATTRS = ('ID', 'Name', 'Alias', 'Parent', 'Target', 'Gap', 'Note')
GFF_ATTRS = dict((a.lower(), a) for a in GFF_ATTRS)

def line_to_dict(line, sep, header_map):
    line = line.rstrip().split(sep)
    assert len(line) == len(header_map), (sep == "\t", line, header_map)
    d = {'attrs':{}}
    for i, value in enumerate(line):
        attrname = header_map[i]
        if attrname in GFF_HEADERS:
            if attrname == 'strand':
                value = '-' if '-' in value else '+'
            d[attrname] = GFF_HEADERS[attrname](value)
        else:
            # it's 'id' try to get 'ID' if it's available in
            # GFF_ATTRS
            d['attrs'][GFF_ATTRS.get(attrname, attrname)] = value
    return d


def to_gff_lines(txtpath, as_dict=False):

    txt = hasattr(txtpath, 'read') and txtpath or open(txtpath)
    header = txt.readline()
    sep = get_sep(header)
    hmap = get_header_map(header, sep)
    for line in txt:
        if line[0] == '#': continue
        fdict = line_to_dict(line, sep, hmap)
        if as_dict:
            yield fdict
        else:
            yield feature_dict_to_gff(fdict)

def feature_dict_to_gff(fdict):
    d = fdict.copy()
    fmt = "%(seqid)s\t%(source)s\t%(type)s\t%(start)i\t%(end)i\t%(score)s"
    fmt += "\t%(strand)s\t%(phase)s\t%(attrs)s"

    for c in ('type', 'source', 'phase', 'score'):
        if not c in d: d[c] = '.'
    attrs = ";".join("%s=%s" % pair for pair in sorted(d['attrs'].items()) if
                     pair[1])
    d['attrs'] = attrs
    s = fmt % d
    return s

if __name__ == "__main__":
    import optparse
    parser = optparse.OptionParser("""
           %s -o out.gff features.txt\n
       where features.txt must contain a header with (at least)
       columns of start, stop/end, seqid/chr/chromosome
       and optionally: name/id, strand, score, phase.
       any additional columns will be used as gff attributes
       """ % sys.argv[0])
    parser.add_option('-o', dest='out_file',
                      help="""path of output gff file if not sent in,\n\
                      the output will go to stdout""")
                       
    (options, args) = parser.parse_args()

    if not args or len(args) > 1:
        print "running tests..."
        import doctest
        doctest.testmod()
        print "\n\n"
        sys.exit(parser.print_help())

    o = options.out_file
    if o:
        o = open(o, 'wb')
    else: o = sys.stdout
    print >>o, GFF_VERSION_HEADER
    for gff_line in to_gff_lines(args[0]):
        print >>o, gff_line
