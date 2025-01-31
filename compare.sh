#!/usr/bin/env sh

IN="$1"
parts=( part.tsv part.xml )
wholes=( whole.tsv whole.xml )

BATCH=320

echo ./iscan.py 12 "$BATCH" "Pieces" "$IN" "${parts[0]}"
     ./iscan.py 12 "$BATCH" "Pieces" "$IN" "${parts[0]}"

sort -n "${parts[0]}" > "${parts[0]}.sorted"

echo ./iscan.sh "$IN" "${wholes[0]}" "${wholes[1]}"
     ./iscan.sh "$IN" "${wholes[0]}" "${wholes[1]}"

sort -n "${wholes[0]}" > "${wholes[0]}.sorted"
