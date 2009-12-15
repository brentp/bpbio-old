given a full genome self-blast and a fasta file. run
::

  python mask_genome.py org_org.blast org.fasta orgname
  python mask_genome.py rice_rice.blast rice.fasta rice

and a new file: rice.masked.fasta will be created in which
any basepair occuring more then 50 times in the blast output
will be masked with 'X'

blast
=====
typical parameters used for the blast are
::

  $ bblast.py -i $FASTA -d $FASTA -a 4 -e 0.001 -m 8 -o ${ORG}_${ORG}.blast

