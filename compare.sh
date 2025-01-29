#!/usr/bin/env sh

IN="$1"
parts=( part.tsv part.xml )
wholes=( whole.tsv whole.xml )

BATCH=320

echo ./iscan.py 12 /tmp "$BATCH" "$IN" "${parts[0]}" "${parts[1]}"
./iscan.py 12 /tmp "$BATCH" "$IN" "${parts[0]}" "${parts[1]}"

echo ./iscan.sh "$IN" "${wholes[0]}" "${wholes[1]}"
./iscan.sh "$IN" "${wholes[0]}" "${wholes[1]}"

# diff -s "${parts[1]}" "${wholes[1]}"
# diff -s "${parts[2]}" "${wholes[2]}"
