#!/bin/bash

file="test_103by300_mixed_matrix.afm"

python test/python/ampute.py $file 0.1 tmp/foo_1.afm
python test/python/ampute.py $file 0.2 tmp/foo_2.afm
python test/python/ampute.py $file 0.3 tmp/foo_3.afm

bin/rf-ace --trainData $file         -i 0 -S 1 -n 100 -m 30 -V tmp/forest_0b.sf -N
bin/rf-ace --trainData tmp/foo_1.afm -i 0 -S 2 -n 100 -m 30 -V tmp/forest_1b.sf -N
bin/rf-ace --trainData tmp/foo_2.afm -i 0 -S 3 -n 100 -m 30 -V tmp/forest_2b.sf -N
bin/rf-ace --trainData tmp/foo_3.afm -i 0 -S 4 -n 100 -m 30 -V tmp/forest_3b.sf -N
bin/rf-ace --trainData $file         -i 0 -S 5 -n 100 -m 30 -V tmp/forest_0t.sf
bin/rf-ace --trainData tmp/foo_1.afm -i 0 -S 6 -n 100 -m 30 -V tmp/forest_1t.sf
bin/rf-ace --trainData tmp/foo_2.afm -i 0 -S 7 -n 100 -m 30 -V tmp/forest_2t.sf
bin/rf-ace --trainData tmp/foo_3.afm -i 0 -S 8 -n 100 -m 30 -V tmp/forest_3t.sf


grep "NNODES=" tmp/forest_0b.sf | cut -d',' -f2 | cut -d'=' -f2 | awk 'BEGIN{ORS="\t"}1;END{print "\n"}' >  tmp/treesizes.tsv
grep "NNODES=" tmp/forest_1b.sf | cut -d',' -f2 | cut -d'=' -f2 | awk 'BEGIN{ORS="\t"}1;END{print "\n"}' >> tmp/treesizes.tsv
grep "NNODES=" tmp/forest_2b.sf | cut -d',' -f2 | cut -d'=' -f2 | awk 'BEGIN{ORS="\t"}1;END{print "\n"}' >> tmp/treesizes.tsv
grep "NNODES=" tmp/forest_3b.sf | cut -d',' -f2 | cut -d'=' -f2 | awk 'BEGIN{ORS="\t"}1;END{print "\n"}' >> tmp/treesizes.tsv
grep "NNODES=" tmp/forest_0t.sf | cut -d',' -f2 | cut -d'=' -f2 | awk 'BEGIN{ORS="\t"}1;END{print "\n"}' >> tmp/treesizes.tsv
grep "NNODES=" tmp/forest_1t.sf | cut -d',' -f2 | cut -d'=' -f2 | awk 'BEGIN{ORS="\t"}1;END{print "\n"}' >> tmp/treesizes.tsv
grep "NNODES=" tmp/forest_2t.sf | cut -d',' -f2 | cut -d'=' -f2 | awk 'BEGIN{ORS="\t"}1;END{print "\n"}' >> tmp/treesizes.tsv
grep "NNODES=" tmp/forest_3t.sf | cut -d',' -f2 | cut -d'=' -f2 | awk 'BEGIN{ORS="\t"}1;END{print "\n"}' >> tmp/treesizes.tsv

 