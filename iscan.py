#!/usr/bin/env python

import datetime
import os
import shlex
import subprocess as sp
import sys
import xml.etree.ElementTree as ET
from multiprocessing import Process
from pathlib import Path
from xml.etree.ElementTree import ElementTree

from Bio import SeqIO

CPUS = int(sys.argv[1])  # os.cpu_count()
BATCH_SIZE = int(sys.argv[2])
PIECES_DIR = Path(sys.argv[3])

INPUT_FAA = Path(sys.argv[4])
OUT_TSV = Path(sys.argv[5])

TMP_FAA = PIECES_DIR / "tmp.faa"

LOG = None
VERBOSE = True
RM_PIECES = False
RM_TMP = False

ENCODING = "UTF-8"


def run_log(cmd, input_stream, log=None, verbose=False):

    def write(msg, handler):
        handler.write(msg)

    if log:
        hlog = open(log, "a", encoding=ENCODING)

    start = datetime.datetime.today()

    msg = f"Starting at {start}\n" + f"Running:\n{shlex.join(cmd)}\n"
    if verbose:
        write(msg, sys.stderr)
    if log:
        write(msg, hlog)

    with open(input_stream, "r", encoding=ENCODING) as input_handler:
        completed = sp.run(
            cmd,
            input=input_handler.read(),
            check=True,
            capture_output=True,
            encoding=ENCODING,
        )

    finish = datetime.datetime.today()

    msg = (
        f"Finished at {finish}\n"
        + f"Time taken was {finish-start}\n"
        + f"STDOUT:\n"
        + f"{completed.stdout}\n"
        + f"STDERR:\n"
        + f"{completed.stderr}\n"
    )
    if verbose:
        write(msg, sys.stderr)

    if log:
        write(msg, hlog)
        hlog.close()


def merge_ETs(*trees: tuple[ElementTree, ...]) -> ElementTree:

    main_tree = trees[0]
    base_root = main_tree.getroot()

    for tree in trees[1:]:
        base_root.extend(tree.getroot())

    return main_tree


def create_fifo(path):
    path = Path(path)
    if not path.exists():
        os.mkfifo(path)
    if not path.is_fifo():
        raise FileExistsError(f"{path} already exists and it is not a FIFO.")
    return path


def xml_cmd_gen(output):
    return [
        "interproscan.sh",
        "--formats",
        "XML",
        "--input",
        "-",
        "--outfile",
        f"{output}",
        "--cpu",
        f"{CPUS}",
        "--goterms",
    ]


def tsv_cmd_gen(output):
    return [
        "interproscan.sh",
        "--mode",
        "convert",
        "--formats",
        "TSV",
        "--input",
        f"-",
        "--outfile",
        f"{output}",
        "--goterms",
        "--enable-tsv-residue-annot",
    ]


def SeqIO_to_xml(seqs, tmp_faa, tmp_xml, log=None, verbose=False) -> None:

    def faa_writer():
        with open(tmp_faa, "w", encoding=ENCODING) as hfaa:
            SeqIO.write(seqs, hfaa, "fasta")

    p = Process(target=faa_writer)
    p.start()

    xml_cmd = xml_cmd_gen(tmp_xml)
    run_log(xml_cmd, tmp_faa, log, verbose)

    p.join()


def is_not_empty(path):
    try:
        return os.path.getsize(path) > 0
    except FileNotFoundError:
        return False


def xmls_to_tsv(xmls, cpus=CPUS):
    pass


if __name__ == "__main__":

    assert INPUT_FAA.exists(), f"File Error: {INPUT_FAA} does NOT exists."
    PIECES_DIR.mkdir(exists_ok=True, parents=True)

    in_faa = SeqIO.parse(INPUT_FAA, format="fasta")
    tmp_faa = create_fifo(TMP_FAA)
    seq_batch = []
    xmls = []

    ixml = 0
    for idx, seq in enumerate(in_faa):
        seq_batch.append(seq)

        if (idx + 1) % BATCH_SIZE == 0:
            ixml += i
            out_xml = PIECES_DIR / f"{ixml}.xml"

            SeqIO_to_XML(seq_batch, TMP_FAA, out_xml, LOG, VERBOSE)
            assert is_not_empty(out_xml)

            xmls.append(out_xml)
            seq_batch = []

    if len(seq_batch) > 0:
        ixml += i
        out_xml = PIECES_DIR / f"{ixml}.xml"

        SeqIO_to_XML(seq_batch, TMP_FAA, out_xml, LOG, VERBOSE)
        assert is_not_empty(out_xml)

        xmls.append(out_xml)
        del seq_batch

    tsv_cmd = tsv_cmd_gen(OUT_TSV)
    run_log(tsv_cmd, OUT_XML, LOG, VERBOSE)

    if RM_TMP:
        TMP_FAA.unlink()

    if RM_PIECES:
        shutil.rmtree(PIECES_DIR)
