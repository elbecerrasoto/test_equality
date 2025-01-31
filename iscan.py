#!/usr/bin/env python

import datetime
import os
import shlex
import shutil
import subprocess as sp
import sys
from multiprocessing import Pool, Process
from pathlib import Path

from Bio import SeqIO

DEBUG = sys.argv[1] == "True"
BATCH_SIZE = int(sys.argv[2])

INPUT_FAA = Path(sys.argv[3])
OUT_TSV = Path(sys.argv[4])

CPUS = os.cpu_count()
PIECES_DIR = Path("pieces")

DEBUG_PATH = PIECES_DIR / Path("debug.sh")

RM_FAA = True
RM_DEBUG = True
RM_PIECES = False

LOG = None
VERBOSE = True

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


if DEBUG:

    DEBUG_SCRIPT = """
    #!/usr/bin/env sh
    OUT="$1"
    cat /dev/stdin > "$OUT"
    """

    with open(DEBUG_PATH, "w") as debugh:
        debugh.write(DEBUG_SCRIPT)

    def xml_cmd_gen(output):
        return ["bash", str(DEBUG_PATH), str(output)]

    def tsv_cmd_gen(output):
        return ["bash", str(DEBUG_PATH), str(output)]

else:

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
            "-",
            "--outfile",
            f"{output}",
            "--goterms",
            "--enable-tsv-residue-annot",
        ]


def SeqIO_to_xml(seqs, tmp_faa, tmp_xml, log=None, verbose=False) -> None:

    with open(tmp_faa, "w", encoding=ENCODING) as hfaa:
        SeqIO.write(seqs, hfaa, "fasta")

    xml_cmd = xml_cmd_gen(tmp_xml)
    run_log(xml_cmd, tmp_faa, log, verbose)

    assert tmp_xml.exists(), f"Non-existent {tmp_xml}"
    assert tmp_faa.exists(), f"Non-existent {tmp_faa}"

    return tmp_xml


def worker_xmls_to_tsvs(input_stream, output_stream):
    cmd = tsv_cmd_gen(output_stream)
    run_log(cmd, input_stream, log=LOG, verbose=VERBOSE)
    return Path(output_stream)


def xmls_to_tsvs(*xmls, cpus=CPUS):
    xmls = [Path(xml) for xml in xmls]
    tsvs = [f"{xml.parent / xml.stem}.tsv" for xml in xmls]
    input_output_pairs = zip(xmls, tsvs)

    with Pool(cpus) as pool:
        results = pool.starmap(worker_xmls_to_tsvs, input_output_pairs)

    for result in results:
        assert result.exists(), f"Non-existent {result}"
    return results


def cat(*files, output=OUT_TSV) -> None:

    with open(output, "wb") as hout:
        for file in files:
            with open(file, "rb") as hfile:
                shutil.copyfileobj(hfile, hout)


if __name__ == "__main__":

    assert INPUT_FAA.exists(), f"File Error: {INPUT_FAA} does NOT exists."
    PIECES_DIR.mkdir(exist_ok=True, parents=True)

    in_faa = SeqIO.parse(INPUT_FAA, format="fasta")
    seq_batch = []
    xmls = []

    ibatch = 0
    for idx, seq in enumerate(in_faa):
        seq_batch.append(seq)

        if (idx + 1) % BATCH_SIZE == 0:
            ibatch += 1

            out_faa = PIECES_DIR / f"{ibatch}.faa"
            out_xml = PIECES_DIR / f"{ibatch}.xml"

            pxml = SeqIO_to_xml(seq_batch, out_faa, out_xml, LOG, VERBOSE)
            xmls.append(pxml)
            seq_batch = []

    if len(seq_batch) > 0:
        ibatch += 1

        out_faa = PIECES_DIR / f"{ibatch}.faa"
        out_xml = PIECES_DIR / f"{ibatch}.xml"

        pxml = SeqIO_to_xml(seq_batch, out_faa, out_xml, LOG, VERBOSE)

        xmls.append(pxml)
        del seq_batch

    tsvs = xmls_to_tsvs(*xmls)
    cat(*tsvs)

    if DEBUG and RM_DEBUG:
        DEBUG_PATH.unlink()

    if RM_PIECES:
        shutil.rmtree(PIECES_DIR)
