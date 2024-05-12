from dataclasses import dataclass
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
    # Dictionary to store fasta sequences
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
    # Dictionaries to store neighboring genes and features
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

            # Check if the line contains "CDS"
            cds_in_line = [wrd for wrd in line if "CDS" in wrd]
            if cds_in_line:
                # Make new variable line_gene to not overwrite line
                line_gene = line

                # Check if the line contains "/gene="
                matching_gene = [wrd for wrd in line_gene if "/gene=" in wrd]

                # While we do not find line with "/gene=" loop through lines.
                while not matching_gene:
                    line_gene = next(gbk_file).strip().split()
                    matching_gene = [wrd for wrd in line_gene if "/gene=" in wrd]

                    # Find "/locus_tag" if there is no "/gene" in the "CDS"
                    # Exit loop and enter elif condition
                    locus_in_line_gene = [
                        wrd for wrd in line_gene if "/locus_tag=" in wrd
                    ]
                    if locus_in_line_gene:
                        break
                if matching_gene:
                    # Save gene name
                    current_gene = [
                        chr.split("=")[1].strip('""') for chr in matching_gene
                    ]

                    # Make new variables, so our while loop would not overwrite variable line_raw
                    line_seq_strip = line_raw.strip()
                    line_seq = line_raw.strip().split()

                    # Find the translation
                    while not [wrd for wrd in line_seq if "/translation=" in wrd]:
                        line_seq_raw = next(gbk_file)
                        line_seq_strip = line_seq_raw.strip()
                        line_seq = line_seq_raw.strip().split()
                    else:
                        if line_seq_strip.endswith('"'):
                            matching_seq = [
                                wrd for wrd in line_seq if "/translation=" in wrd
                            ]
                            current_seq = [
                                chr.split("=")[1].strip('""') for chr in matching_seq
                            ]

                        elif not line_seq_strip.endswith('"'):
                            while not line_seq_strip.endswith('"'):
                                line_seq_strip = next(gbk_file).strip()
                                current_seq += line_seq_strip

                        current_gene = "".join(map(str, current_gene))
                        current_seq = "".join(map(str, current_seq))
                        features[current_gene] = current_seq

                # Condition for "/locus_tag="
                elif locus_in_line_gene:
                    current_gene = [
                        chr.split("=")[1].strip('""') for chr in locus_in_line_gene
                    ]

                    # Make new variables, so our while loop would not overwrite variable line_raw
                    line_seq_raw = line_raw.strip()
                    line_seq = line_raw.strip().split()

                    # Find the translation
                    while not [wrd for wrd in line_seq if "/translation=" in wrd]:
                        line_seq_raw = next(gbk_file)
                        line_seq_strip = line_seq_raw.strip()
                        line_seq = line_seq_raw.strip().split()
                    else:
                        if line_seq_strip.strip().endswith('"'):
                            matching_seq = [
                                wrd for wrd in line_seq if "/translation=" in wrd
                            ]
                            current_seq = [
                                chr.split("=")[1].strip('""') for chr in matching_seq
                            ]

                        elif not line_seq_strip.strip().endswith('"'):
                            while not line_seq_strip.strip().endswith('"'):
                                line_seq_strip = next(gbk_file).strip()
                                current_seq += line_seq_strip

                        current_gene = "".join(map(str, current_gene))
                        current_seq = "".join(map(str, current_seq))
                        features[current_gene] = current_seq

    # Extract neighboring genes
    for gene in genes:
        for i, k in enumerate(features.keys()):
            if gene in k:
                for j in range(n_before):
                    if i - (j + 1) >= 0:
                        gene_key, seq_val = list(features.items())[i - (j + 1)]
                        neighbours_before[gene_key] = seq_val
                for j in range(n_after):
                    if i + (j + 1) < len(features):  # Fix off-by-one error
                        gene_key, seq_val = list(features.items())[i + (j + 1)]
                        neighbours_after[gene_key] = seq_val

    # Determine the output FASTA file path
    if output_fasta is None:
        fasta_dir = os.path.dirname(input_gbk)
        output_fasta = (
            Path(input_gbk).with_suffix("").stem + "_flanking_genes" + ".fasta"
        )
    fasta_dir = os.path.dirname(input_gbk)
    fasta_path = os.path.join(fasta_dir, output_fasta)

    # Write the extracted sequences to the output FASTA file
    with open(fasta_path, "w") as fasta_file:
        for gene, sequence in neighbours_before.items():
            fasta_file.write(f">{gene}\n{sequence}\n")
        for gene, sequence in neighbours_after.items():
            fasta_file.write(f">{gene}\n{sequence}\n")


@dataclass
class FastaRecord:
    id: str  # Unique identifier for the sequence
    seq: str  # The actual DNA or protein sequence
    description: str  # Description associated with the sequence (optional)

    # Define a custom string representation for the FastaRecord object
    def __repr__(self):
        return (
            f"ID: {self.id},\n Description: {self.description},\n Sequence:{self.seq}\n"
        )


class OpenFasta:
    def __init__(self, file_path: str):
        """
        Initialize an OpenFasta object with the file path.

        Args:
            file_path (str): Path to the FASTA file.
        """
        self.file_path = file_path
        self.current_record: str = None  # Stores the current FASTA header line
        self.header = None  # Flag to indicate if currently processing header
        self.stop = False  # Flag to indicate end of file

    # Context manager support for using OpenFasta with `with` statement
    def __enter__(self):
        self.handler = open(self.file_path)
        return self

    # Define the iterator protocol for OpenFasta
    def __iter__(self):
        return self

    # Generator function to yield FastaRecord objects one by one
    def __next__(self):
        seq = ""
        header = self.header
        stop = self.stop
        current_record = self.current_record

        if stop:
            raise StopIteration("End of file")

        # Process existing record if available
        elif current_record is not None:
            line = current_record
            id, description = line.split(" ", 1)  # Split header line by first space
            header = True

        # Read and process a new record if no current record
        elif header is None:
            line = self.handler.readline().strip()
            id, description = line.split(" ", 1)  # Split header line by first space
            header = True

        # Read sequence lines until next header or end of file
        while header and not stop:
            line = self.handler.readline().strip()
            if line == "":
                self.stop = True
                return FastaRecord(id, seq, description)
            elif line[0] == ">":
                self.current_record = line
                return FastaRecord(id, seq, description)
            else:
                seq += line

    # Provide convenience methods to read single or all records
    def read_record(self):
        return self.__next__()

    def read_records(self):
        return list(self.__iter__())

    # Close the file when exiting the context manager
    def __exit__(self, exc_type, exc_val, exc_tb):
        self.handler.close()
