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

DEBUG = bool(sys.argv[1])
BATCH_SIZE = int(sys.argv[2])

INPUT_FAA = Path(sys.argv[3])
OUT_TSV = Path(sys.argv[4])

CPUS = os.cpu_count()  # int(sys.argv[2])
PIECES_DIR = Path("pieces")

TMP_FAA = PIECES_DIR / "tmp.faa"

LOG = None
VERBOSE = True
RM_PIECES = False
RM_TMP = True

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


def create_fifo(path):
    path = Path(path)
    if not path.exists():
        os.mkfifo(path)
    if not path.is_fifo():
        raise FileExistsError(f"{path} already exists and it is not a FIFO.")
    return path


if DEBUG:

    DEBUG_PATH = Path("debug.sh")
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
    tmp_faa = create_fifo(TMP_FAA)
    seq_batch = []
    xmls = []

    ixml = 0
    for idx, seq in enumerate(in_faa):
        seq_batch.append(seq)

        if (idx + 1) % BATCH_SIZE == 0:
            ixml += 1
            out_xml = PIECES_DIR / f"{ixml}.xml"

            SeqIO_to_xml(seq_batch, TMP_FAA, out_xml, LOG, VERBOSE)
            assert out_xml.exists(), f"Non-existent {out_xml}"

            xmls.append(out_xml)
            seq_batch = []

    if len(seq_batch) > 0:
        ixml += 1
        out_xml = PIECES_DIR / f"{ixml}.xml"

        SeqIO_to_xml(seq_batch, TMP_FAA, out_xml, LOG, VERBOSE)
        assert out_xml.exists(), f"Non-existent {out_xml}"

        xmls.append(out_xml)
        del seq_batch

    tsvs = xmls_to_tsvs(*xmls)
    cat(*tsvs)

    if RM_TMP:
        TMP_FAA.unlink()

    if RM_PIECES:
        shutil.rmtree(PIECES_DIR)
