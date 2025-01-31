#!/usr/bin/env sh

fasta_sort < out.faa > out.s.faa
fasta_sort < all.faa > all.s.faa

diff -s out.s.faa all.s.faa
