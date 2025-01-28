#!/usr/bin/env sh

IN="$1"
parts=( part.tsv part.xml )
wholes=( whole.tsv whole.xml )

BATCH=320

./iscan.py 12 /tmp "$BATCH" "$IN" "${parts[1]}" "${parts[2]}"
./iscan.sh "$IN" "${wholes[1]}" "${wholes[2]}"

# diff -s "${parts[1]}" "${wholes[1]}"
# diff -s "${parts[2]}" "${wholes[2]}"
