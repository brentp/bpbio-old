==============
Genome Masking
==============

given a full genome self-blast and a fasta file, create
a new masked genome where any basepair occuring more than
`cutoff` times is masked. the new file will be written
in the same directory as the original fasta with a name like:
    rice.masked.50.fasta

blast
=====
typical parameters used for the blast are
::

  $ bblast.py -i $FASTA -d $FASTA -a 4 -e 0.001 -m 8 -o ${ORG}_${ORG}.blast

shell
=====
an entire mask. start from just the genomic fasta can be run with
::

    ORG=thaliana_v7
    FA=${ORG}.fasta
    BL=${ORG}_${ORG}.blast
    bblast.py -p blastn -i $FA -d $FA -m 8 -a 8 -o $BL
    python mask_genome.py -b $BL -f $FA -o $ORG -c 50

after running the self-self blast on the thaliana_v7.fasta,
this will create a new fasta file "thaliana_v7.masked.50.fasta"
with all basepairs covered by more than 50 blast hits masked to 'X'


todo
====
1. allow soft-masking.
