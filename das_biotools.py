from abc import ABC, abstractmethod
from dataclasses import dataclass

import os
import sys
import datetime
import requests
import io
from bs4 import BeautifulSoup
from typing import List
from dotenv import load_dotenv

load_dotenv()
TG_API_TOKEN = os.getenv("TG_API_TOKEN")


# HW 14
class BiologicalSequence(ABC):
    @abstractmethod
    def __init__(self, sequence):
        pass

    @abstractmethod
    def __len__(self):
        pass

    @abstractmethod
    def __getitem__(self, slc):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __repr__(self):
        pass

    @abstractmethod
    def is_valid_alphabet(self):
        pass


class NucleicAcid(BiologicalSequence):
    def __init__(self, sequence):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, slc):
        return self.sequence[slc]

    def __str__(self):
        return str(self.sequence)

    def __repr__(self):
        return self.sequence

    def is_valid_alphabet(self):
        alphabet = type(self).ALPHABET
        if set(self.sequence).issubset(alphabet):
            return True
        else:
            return False

    def complement(self):
        if type(self) == NucleicAcid:
            raise NotImplementedError("Cannot complement NucleicAcid instance")

        map_dict = type(self).MAP
        comp_seq = "".join([map_dict[base] for base in self.sequence])

        return type(self)(comp_seq)

    def gc_content(self):
        if type(self) == NucleicAcid:
            raise NotImplementedError("Cannot gc_content NucleicAcid instance")
        gc_count = sum([1 for base in self.sequence if base in ["C", "G"]])
        gc_content = (gc_count / len(self)) * 100

        return gc_content


class DNASequence(NucleicAcid):
    ALPHABET = set("ATGC")
    MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}

    def transcribe(self):
        transcribed = self.sequence.replace("T", "U")

        return RNASequence(transcribed)


class RNASequence(NucleicAcid):
    ALPHABET = set("AUGC")
    MAP = {"A": "U", "U": "A", "C": "G", "G": "C"}


class AminoAcidSequence(BiologicalSequence):
    ALPHABET = set("ACDEFGHIKLMNPQRSTVWY")

    def __init__(self, sequence):
        self.sequence = sequence

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, slc):
        return self.sequence[slc]

    def __str__(self):
        return str(self.sequence)

    def __repr__(self):
        return self.sequence

    def is_valid_alphabet(self):
        alphabet = type(self).ALPHABET
        if set(self.sequence).issubset(alphabet):
            return True
        else:
            return False

    amino_acid_frequency = {}

    def calculate_aa_freq(self):
        """
        Calculates the frequency of each amino acid in a protein sequence or sequences.

        :param sequences: protein sequence or sequences
        :type sequences: str or list of str
        :return: dictionary with the frequency of each amino acid
        :rtype: dict
        """

        # Creating a dictionary with aminoacid frequencies:
        amino_acid_frequency = {}

        for amino_acid in self.sequence:
            # If the aminoacid has been already in:
            if amino_acid in amino_acid_frequency:
                amino_acid_frequency[amino_acid] += 1
            # If the aminoacid hasn't been already in:
            else:
                amino_acid_frequency[amino_acid] = 1

        return amino_acid_frequency


# HW 17


def telegram_logger(chat_id):
    def decorator(func):
        def wrapper(*args, **kwargs):
            default_stdout = sys.stdout
            default_stderr = sys.stderr
            # Redirect stdout and stderr to StringIO
            stdout_stderr_file = io.StringIO()
            sys.stdout = sys.stderr = stdout_stderr_file

            start_time = datetime.datetime.now()
            try:
                func(*args, **kwargs)
                end_time = datetime.datetime.now()
                duration = end_time - start_time
                if duration < datetime.timedelta(days=1):
                    duration_str = str(duration)
                else:
                    days = duration.days
                    hours, remainder = divmod(duration.seconds, 3600)
                    minutes, seconds = divmod(remainder, 60)
                    duration_str = f"{days} days, {hours:02}:{minutes:02}:{seconds:02}"
                message = f"Function `{func.__name__}` executed successfully in `{duration_str}`."
                send_telegram_message(
                    chat_id, message, stdout_stderr_file.getvalue(), good_function=True
                )
            except Exception as e:
                end_time = datetime.datetime.now()
                duration = end_time - start_time
                message = f"Function `{func.__name__}` execution failed with error: \n `{type(e).__name__}: {e}`"
                send_telegram_message(
                    chat_id, message, stdout_stderr_file.getvalue(), good_function=False
                )
                raise
            finally:
                # Reset stdout and stderr
                sys.stdout = default_stdout
                sys.stderr = default_stderr

        return wrapper

    return decorator


def send_telegram_message(chat_id, message, stdout="", good_function=None):
    url = f"https://api.telegram.org/bot{TG_API_TOKEN}/sendDocument"

    params = {"chat_id": chat_id, "caption": message, "parse_mode": "Markdown"}

    files = {}
    if stdout:
        if good_function:
            files["document"] = ("good_function.log", stdout)
        else:
            files["document"] = ("bad_function.log", stdout)
    try:
        response = requests.post(url, params=params, files=files)
        response.raise_for_status()
        print("Message sent successfully to Telegram!")
    except requests.exceptions.RequestException as e:
        print(f"Error sending message to Telegram: {e}")


@dataclass
class GenscanOutput:
    """
    Represents the output of the GENSCAN prediction.

    Attributes:
        status (int): The HTTP status code of the response.
        cds_list (List[dict]): List of predicted peptides.
        intron_list (List[dict]): List of predicted introns.
        exon_list (List[dict]): List of predicted exons.
    """

    status: int
    cds_list: List[dict]
    intron_list: List[dict]
    exon_list: List[dict]


def run_genscan(
    sequence=None,
    sequence_file=None,
    organism="Vertebrate",
    exon_cutoff=1.00,
    sequence_name="",
):
    """
    Run GENSCAN prediction.

    Args:
        sequence (str, optional): The input DNA sequence.
        sequence_file (str, optional): Path to the file containing the DNA sequence.
        organism (str, optional): The organism type for prediction. Defaults to "Vertebrate".
        exon_cutoff (float, optional): Exon cutoff value. Defaults to 1.00.
        sequence_name (str, optional): Name of the sequence. Defaults to "".

    Returns:
        GenscanOutput: An instance of GenscanOutput containing prediction results.
    """

    # Read sequence from file if provided
    if sequence_file:
        with open(sequence_file, "r") as file:
            sequence = file.read().strip()

    # Prepare parameters for the POST request
    params = {
        "-o": organism,
        "-s": sequence,
        "-n": sequence_name,
        "-e": exon_cutoff,
        "-p": "Predicted peptides only",
    }

    # Send POST request to GENSCAN service
    url = "http://argonaute.mit.edu/cgi-bin/genscanw_py.cgi"
    response = requests.post(url, data=params)
    status = response.status_code

    # Check if connection was successful
    if status != 200:
        raise ValueError("Failed to connect to GENSCAN service")

    # Parse the HTML response
    soup = BeautifulSoup(response.text, "html.parser")
    pre_tag = soup.find("pre")
    pre_text = pre_tag.get_text()
    lines = pre_text.split("\n")

    # Extract predicted exons
    exon_list = extract_exons(lines)

    # Extract predicted introns
    intron_list = extract_introns(exon_list)

    # Extract predicted peptides
    cds_list = extract_predicted_peptides(lines)

    return GenscanOutput(status, cds_list, intron_list, exon_list)


def extract_exons(lines):
    """
    Extract predicted exons from the lines of text.

    Args:
        lines (List[str]): Lines of text from the GENSCAN prediction.

    Returns:
        List[dict]: List of predicted exons.
    """
    exon_list = []
    exon_start_section = lines.index(
        "Gn.Ex Type S .Begin ...End .Len Fr Ph I/Ac Do/T CodRg P.... Tscr.."
    )
    exon_end_section = lines.index("Suboptimal exons with probability > 1.000")
    exon_lines = lines[exon_start_section:exon_end_section]

    exon_section = None
    for line in exon_lines:
        if line.startswith(" 1"):
            exon_section = True
        if exon_section and line:
            parts = line.split()
            number = parts[0]
            type = parts[1]
            start = int(parts[3])
            end = int(parts[4])
            exon_info = {"number": number, "type": type, "start": start, "end": end}
            exon_list.append(exon_info)
    return exon_list


def extract_introns(exon_list):
    """
    Extract predicted introns from the list of predicted exons.

    Args:
        exon_list (List[dict]): List of predicted exons.

    Returns:
        List[dict]: List of predicted introns.
    """
    intron_list = []
    for i in range(len(exon_list) - 1):
        current_exon = exon_list[i]
        next_exon = exon_list[i + 1]

        intron_start = None
        intron_end = None

        if current_exon["end"] < current_exon["start"]:
            if next_exon["end"] < next_exon["start"]:
                intron_start = current_exon["start"] + 1
                intron_end = next_exon["end"] - 1
            else:
                intron_start = current_exon["start"] + 1
                intron_end = next_exon["start"] - 1
        else:
            if next_exon["end"] < next_exon["start"]:
                intron_start = current_exon["end"] + 1
                intron_end = next_exon["end"] - 1
            else:
                intron_start = current_exon["end"] + 1
                intron_end = next_exon["start"] - 1

        intron_list.append(
            {
                "number": f"{current_exon['number']}-{next_exon['number']}",
                "type": "Intron",
                "start": intron_start,
                "end": intron_end,
            }
        )
    return intron_list


def extract_predicted_peptides(lines):
    """
    Extract predicted peptides from the lines of text.

    Args:
        lines (List[str]): Lines of text from the GENSCAN prediction.

    Returns:
        List[dict]: List of predicted peptides.
    """
    cds_list = []
    # Find start index of predicted peptides section
    cds_start_section = lines.index("Predicted peptide sequence(s):")
    cds_lines = lines[cds_start_section:]

    cds_section = None
    cds_seq = []
    for line in cds_lines:
        if line.startswith(">"):
            cds_section = True
            line_split = line.split("|")
            peptide_info = line_split[1:]
            continue

        if cds_section and line:
            cds_seq.append(line)

    cds_list.append(
        {
            "peptide": peptide_info,
            "sequence": "".join(cds_seq),
        }
    )
    return cds_list
