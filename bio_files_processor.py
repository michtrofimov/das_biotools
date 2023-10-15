import os
from pathlib import Path
from typing import List


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None):
    fasta = dict()
    with open(input_fasta, "r") as fasta_file:
        seq = ""
        id = ""
        for i, line in enumerate(fasta_file):
            line = line.strip()
            if line.startswith(">") and i == 0:
                id = line
                fasta[id] = None
            elif line.startswith(">"):
                fasta[id] = seq
                seq = ""
                id = line
            else:
                seq += line

    if output_fasta is None:
        fasta_dir = os.path.dirname(input_fasta)
        output_fasta = Path(input_fasta).with_suffix("").stem + "_converted.fasta"

    fasta_path = os.path.join(fasta_dir, output_fasta)

    with open(fasta_path, "w") as fasta_file:
        for key, value in fasta.items():
            fasta_file.write(key + "\n")
            fasta_file.write(value + "\n")
