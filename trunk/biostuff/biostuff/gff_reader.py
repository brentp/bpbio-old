import mmap
import re
import sys

import logging
logging.basicConfig(level=logging.DEBUG)

def _get_non_comment_line(fh):
    while True:
        line = fh.readline()
        if line and line[0] == '#': continue
        return line

def _get_lines_until_next(fh, parent_ids):
    """
    so this just continues to append lines to lines[]
    until Parent attribute of a line does not match
    any of the values in parent_ids
    any time a line is found with a new ID attribute
    (whose Parent attr matches the current parent_ids list),
    that line's own ID is added to the parent_ids list"""

    lines = [_get_non_comment_line(fh)]
    while True:
        new_parent = False
        for parent_id in parent_ids:
            if 'Parent=' + parent_id in lines[-1]:
                if not lines[-1]: return lines # end of file
                lines.append(_get_non_comment_line(fh))
                if 'ID=' in lines[-1]:
                    new_parent=True
                break
        else:
            break
        if new_parent:
            parent_ids.append(GFFLine(lines[-1]).attribs["ID"])
            parent_ids = list(set(parent_ids))

    return lines

class GFFNode(object):
    def __init__(self, node_list):
        self.start = min(n.start for n in node_list)
        self.stop = max(n.stop for n in node_list)
        self.end = self.stop
        assert "ID" in node_list[0].attribs, (node_list[0], node_list[0].attribs)
        self.parent = node_list[0]
        self.nodes = node_list[1:]

        if self.parent.start != self.start:
            logging.debug(("the start of the parent != the start of the item:" + \
                        ", ".join(map(str, (self, self.start, self.parent.start)))))
        if self.parent.end != self.end:
            logging.debug("the end of the parent != the end of the item:" + \
                    ", ".join(map(str, (self, self.end, self.parent.end))))
    
    def __repr__(self):
        return "GFFNode(%s: %i .. %i, %i sub-nodes)" % \
                   (self.parent.attrs.get('ID', ''), self.start, self.end,
                    len(self.nodes))


    @classmethod
    def yield_nodes(cls, fh):
        close = False
        if not hasattr(fh, 'read'):
            close = True
            fh = open(fh, 'r')

        next_gene_line = _get_non_comment_line(fh)
        while next_gene_line:
            assert "ID=" in next_gene_line,\
                    ("should have an id to be a parent feature",)

            parent = GFFLine(next_gene_line)
            block = _get_lines_until_next(fh, [parent.attribs["ID"]])
            block.insert(0, next_gene_line)
            next_gene_line = block.pop()

            yield GFFNode([GFFLine(l) for l in block])


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

