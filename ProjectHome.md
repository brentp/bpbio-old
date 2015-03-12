bioinformatics stuffs. any feedback or contributions would be **greatly** appreciated i can move code to github/bitbucket to facilitate that if it helps:

**NOTE**: as of 2009/10/13 some modules which previously lived here have been moved to [bitbucket](http://bitbucket.org/brentp/biostuff/)

  * biostuff: fast, memory efficient (cython) reader and objects for GFF and tab-delimited BLAST files: http://code.google.com/p/bpbio/source/browse/#svn/trunk/biostuff/biostuff
```
    from biostuff import BlastLine
    for line in open('some_file.blast'):
        b = BlastLine(line)
        print b.query, b.subject, b.pctid, b.qstart, b.qstop
```

  * gfftools: these mostly require genometools to work, but allow merging and sorting of existing gff files and drawing (using matplotlib) of gff files onto chromosomes.

  * ihash : this is a project to encode an interval to a string suitable for indexing in a b-tree, has the same limitations as geohash (http://en.wikipedia.org/wiki/Geohash), on which it is based.

  * interval\_tree : a pure python interval tree that's extremely simple. also a shedskin version. for real work, use the interval tree in bx-python instead.

  * seqfind : a cython BKTree and damerau-levenshtein distance for indexing words/sequences.

  * scripts/dagchainer : scripts for and modifications to [B. Haas' dagchainer](http://dagchainer.sourceforge.net) for finding syntenic pairs.

  * scripts/localdups : find local (tandem) duplications given a simple input format.

  * scripts/mask\_genome : given a full-genome BLAST output and a fasta file, create a masked version of the fasta based on the number of BLAST hits at each basepair location.

  * gtpym : a module to add some features to the [genometools](http://genometools.org) python library.

  * and more ...