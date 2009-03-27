import mmap
import re
import sys

def get_non_comment_line(fh):
    while True:
        line = fh.readline()
        if line[0] == '#': continue
        return line

def get_lines_until_next(fh, feature):
    lines = [get_non_comment_line(fh)]
    while not feature in lines[-1]:
        if not lines[-1]: return lines # end of file
        lines.append(get_non_comment_line(fh))
    return lines

class GFFNode(object):
    def __init__(self):
        pass

    @classmethod
    def yield_nodes(cls, fh, top_level_feature='gene'):
        gene = top_level_feature
        if not hasattr(fh, 'read'):
            fh = open(fh, 'r')

        next_gene_line = get_non_comment_line(fh)
        while True:
            assert gene in next_gene_line
            block = get_lines_until_next(fh, gene)
            block.insert(0, next_gene_line)
            next_gene_line = block.pop()

            yield [GFFLine(l) for l in block]




class GFFLine(object):
    __slots__ = ('seqid', 'com', 'type', 'start', 'stop', 'end', 'strand',
                    'attrs', 'attribs')
    def __init__(self, sline):
        line = sline.split("\t")
        self.seqid = line[0]
        self.com  = line[1]
        self.type = line[2]
        self.start = int(line[3])
        self.stop = self.end = int(line[4])
        self.strand = line[6] in ('-', '-1') and -1 or 1
        self.attrs = self.attribs = self._parse_attrs(line[8])
    
    def _parse_attrs(self, sattrs):
        attrs = {}
        if "=" in sattrs:
            for pair in sattrs.split(';'):
                if not pair: continue
                pair = pair.split('=')
                attrs[pair[0]] = pair[1]
        if attrs == {}: attrs["ID"] = line[8]
        return attrs

    def __repr__(self):
        return "GFFLine(%s %s:%i .. %i)" % (self.seqid,
                          self.attrs.get("ID", ""), self.start, self.end)

