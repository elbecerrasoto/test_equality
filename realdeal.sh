#!/usr/bin/env sh

./iscan.py False 380 all.faa parts.tsv
./iscan.sh all.faa whole.tsv whole.xml

sort -n parts.tsv > parts.s.tsv
sort -n whole.tsv > whole.s.tsv

diff -s whole.s.tsv parts.s.tsv


