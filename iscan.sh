#!/usr/bin/env sh

IN_FAA="$1"
OUT_TSV="$2"
OUT_XML="$3"

interproscan.sh --formats XML\
    --input "$IN_FAA"\
    --outfile "$OUT_XML"\
    --cpu 12\
    --goterms

interproscan.sh --mode convert\
    --formats TSV\
    --input "$OUT_XML"\
    --outfile "$OUT_TSV"\
    --goterms --enable-tsv-residue-annot
