=======
pyfasta
=======

:Description: pythonic access to fasta sequence files.

:Author: Brent Pedersen (brentp)
:Email: bpederse@gmail.com
:License: MIT


Implementation
==============

Requires Python >= 2.6. Stores a flattened version of the fasta file without 
spaces or headers. And a pickle of the start, stop (for fseek) locations of 
each header in the fasta file for internal use.


Usage
=====


.. sourcecode:: python

    >>> from pyfasta import Fasta

    >>> f = Fasta('tests/data/three_chrs.fasta')
    >>> sorted(f.keys())
    ['chr1', 'chr2', 'chr3']

    >>> f['chr1']
    FastaRecord('tests/data/three_chrs.fasta.flat', 0..80)

    >>> f['chr1'][:10]
    'ACTGACTGAC'


    # the index stores the start and stop of each header from teh fasta file
    >>> f.index
    {'chr3': (160, 3760), 'chr2': (80, 160), 'chr1': (0, 80)}


    # can query by a 'feature' dictionary
    >>> f.sequence({'chr': 'chr1', 'start': 2, 'stop': 9})
    'CTGACTGA'

    # with reverse complement for - strand
    >>> f.sequence({'chr': 'chr1', 'start': 2, 'stop': 9, 'strand': '-'})
    'TCAGTCAG'


    # creates a .flat and a .gdx pickle of the fasta and the index.
    >>> import os
    >>> sorted(os.listdir('tests/data/'))[1:]
    ['three_chrs.fasta', 'three_chrs.fasta.flat', 'three_chrs.fasta.gdx']

    # cleanup (though for real use these will remain for faster access)
    >>> os.unlink('tests/data/three_chrs.fasta.gdx')
    >>> os.unlink('tests/data/three_chrs.fasta.flat')
