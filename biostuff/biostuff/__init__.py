#from blast_line import BlastLine as PyBlastLine
from cblastline import BlastLine
from gff_reader import GFFNode, GFFLine


def main():
    import sys
    progs = ['cblastline']
    if len(sys.argv) == 1 or not sys.argv[1] in progs:
        print """Usage:
    $ biostuff prog [opts]
or 
    $ biostuff prog -h

where current prog's are: %s
        """ % "\n".join(progs)
    else:
        prog = sys.argv.pop(1)
        mod = __import__(prog)
        mod.main()

if __name__ == "__main__":
    main()
