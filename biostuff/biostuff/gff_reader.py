import mmap
import re
import sys

def _get_non_comment_line(fh):
    while True:
        line = fh.readline()
        if line[0] == '#': continue
        return line

def _get_lines_until_next(fh, feature):
    lines = [_get_non_comment_line(fh)]
    while not feature in lines[-1]:
        if not lines[-1]: return lines # end of file
        lines.append(_get_non_comment_line(fh))
    return lines

class GFFNode(object):
    def __init__(self, node_list, top_level_feature='gene'):
        self.start = min(n.start for n in node_list)
        self.stop = max(n.stop for n in node_list)
        self.end = self.stop
        assert node_list[0].type == top_level_feature, node_list[0]
        self.parent = node_list[0]
        self.nodes = node_list[1:]

        assert self.parent.start == self.start, (self.start, self.parent.start)
        assert self.parent.end == self.end, (self.end, self.parent.end)
    
    def __repr__(self):
        return "GFFNode(%s: %i .. %i, %i sub-nodes)" % \
                   (self.parent.attrs.get('ID', ''), self.start, self.end,
                    len(self.nodes))


    @classmethod
    def yield_nodes(cls, fh, top_level_feature='gene'):
        gene = top_level_feature
        close = False
        if not hasattr(fh, 'read'):
            close = True
            fh = open(fh, 'r')

        next_gene_line = _get_non_comment_line(fh)
        while True:
            assert gene in next_gene_line
            block = _get_lines_until_next(fh, gene)
            block.insert(0, next_gene_line)
            next_gene_line = block.pop()

            yield GFFNode([GFFLine(l) for l in block], top_level_feature=gene)


        if close: fh.close()


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

    @classmethod
    def yield_lines(cls, fh):
        close = False
        if not hasattr(fh, 'read'):
            fh = open(fh, 'r')
            close = True
        line = _get_non_comment_line(fh)
        while line:
            yield GFFLine(line)
            line = _get_non_comment_line(fh)
    
        if close: fh.close()

    def __repr__(self):
        return "GFFLine(%s %s:%i .. %i)" % (self.seqid,
                          self.attrs.get("ID", ""), self.start, self.end)

