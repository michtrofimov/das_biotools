import os
from pathlib import Path
from typing import List


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None):
    """
    Convert a multi-line FASTA file to a one-line FASTA file.

    Parameters:
    - input_fasta (str): Path to the input multi-line FASTA file.
    - output_fasta (str): Path to the output one-line FASTA file.
                         If not provided, a default name will be used.

    Returns:
    - FASTA file
    """
    fasta = dict()

    # Read the input multi-line FASTA file
    with open(input_fasta, "r") as fasta_file:
        seq = ""
        id = ""
        for i, line in enumerate(fasta_file):
            line = line.strip()

            # Identify sequence identifiers (lines starting with ">")
            if line.startswith(">") and i == 0:
                id = line
                fasta[id] = None
            elif line.startswith(">"):
                fasta[id] = seq
                seq = ""
                id = line
            else:
                seq += line

    # Determine the output one-line FASTA file path
    if output_fasta is None:
        fasta_dir = os.path.dirname(input_fasta)
        output_fasta = Path(input_fasta).with_suffix("").stem + "_converted.fasta"

    fasta_path = os.path.join(fasta_dir, output_fasta)

    # Write the converted sequences to the output one-line FASTA file
    with open(fasta_path, "w") as fasta_file:
        for key, value in fasta.items():
            fasta_file.write(key + "\n")
            fasta_file.write(value + "\n")


def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: List[str],
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = None,
):
    """
    Extracts gene sequences from a GenBank (gbk) file and creates a FASTA file
    with specified neighboring genes.

    Parameters:
    - input_gbk (str): Path to the input GenBank file.
    - genes (List[str]): List of genes to extract.
    - n_before (int): Number of genes before the specified genes to include.
    - n_after (int): Number of genes after the specified genes to include.
    - output_fasta (str): Path to the output FASTA file.

    Returns:
    - FASTA file
    """
    neighbours_before = dict()
    neighbours_after = dict()
    features = dict()

    # Handle the case when a single gene is provided as a string
    if isinstance(genes, str):
        genes = [genes]

    current_gene = None
    current_seq = ""

    with open(input_gbk, "r") as gbk_file:
        for line_raw in gbk_file:
            line = line_raw.strip().split()

            cds_in_line = any("CDS" in word for word in line)
            if cds_in_line:
                matching_gene = [word for word in line if ("/gene=" in word)]

                while not matching_gene:
                    line = next(gbk_file).strip().split()
                    matching_gene = [word for word in line if ("/gene=" in word)]

                    locus_in_line_gene = [
                        word for word in line if ("/locus_tag=" in word)
                    ]
                    if locus_in_line_gene:
                        break

                if matching_gene:
                    current_gene = "".join(
                        map(
                            str,
                            [word.split("=")[1].strip('""') for word in matching_gene],
                        )
                    )
                    current_seq = ""

                    while not current_seq.strip().endswith(('"')):
                        current_seq += line_raw.strip()
                        line_raw = next(gbk_file).strip()

                    current_seq = current_seq.strip().split("=")[1].strip('""')
                    features[current_gene] = current_seq

                elif locus_in_line_gene:
                    current_gene = "".join(
                        map(
                            str,
                            [
                                word.split("=")[1].strip('""')
                                for word in locus_in_line_gene
                            ],
                        )
                    )
                    current_seq = ""

                    while not current_seq.strip().endswith(('"')):
                        current_seq += line_raw.strip()
                        line_raw = next(gbk_file).strip()

                    current_seq = current_seq.strip().split("=")[1].strip('""')
                    features[current_gene] = current_seq

    # Extract neighboring genes
    for gene in genes:
        for i, (gene_key, seq_val) in enumerate(features.items()):
            if gene in gene_key:
                for j in range(n_before):
                    if i - (j + 1) >= 0:
                        neighbours_before[gene_key] = seq_val
                for j in range(n_after):
                    if i + (j + 1) < len(features):
                        neighbours_after[gene_key] = seq_val

    # Determine the output FASTA file path
    if output_fasta is None:
        fasta_dir = os.path.dirname(input_gbk)
        output_fasta = (
            Path(input_gbk).with_suffix("").stem + "flanking_genes" + ".fasta"
        )
    fasta_path = os.path.join(fasta_dir, output_fasta)

    # Write the extracted sequences to the output FASTA file
    with open(fasta_path, "w") as fasta_file:
        for gene, sequence in neighbours_before.items():
            fasta_file.write(f">{gene}\n{sequence}\n")
        for gene, sequence in neighbours_after.items():
            fasta_file.write(f">{gene}\n{sequence}\n")


# Example usage
select_genes_from_gbk_to_fasta("tests/example_gbk.gbk", "flu_2", n_before=2)
