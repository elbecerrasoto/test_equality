#!/usr/bin/env sh

./iscan.py False 380 all.faa parts.tsv
./iscan.sh all.faa whole.tsv whole.xml

fasta_sort < pieces/1.faa >| pieces/1.s.faa
fasta_sort < all.faa >| all.s.faa
diff -s pieces/1.s.faa all.s.faa


sort -n parts.tsv > parts.s.tsv
sort -n whole.tsv > whole.s.tsv
diff -s whole.s.tsv parts.s.tsv
