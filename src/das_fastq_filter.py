import os
from typing import Dict, Tuple


def read_fastq_file(input_path: str) -> Dict[str, Tuple[str, str]]:
    """
    Read a FASTQ file and return a dictionary where keys are sequence identifiers
    and values are tuples containing sequence and quality scores.

    Parameters:
    - input_path (str): Path to the input FASTQ file.

    Returns:
    - dict: A dictionary with sequence identifiers as keys and tuples (sequence, quality scores) as values.
    """

    fastq = dict()
    with open(input_path, "r") as fastq_file:
        lines = list()
        for line in fastq_file:
            lines.append(line.strip())
            if len(lines) == 4:
                id = lines[0]
                (read, qual) = lines[1::2]
                fastq[id] = (read, qual)
                lines = list()

    return fastq


def calculate_gc_content(seq: str) -> float:
    """
    Calculate a GC content percentage of a DNA sequence.

    Parameters:
    - seq (str): DNA sequence

    Returns:
    - float: GC content percentage
    """

    gc_content = 0

    for nt in seq:
        if (nt == "G") or (nt == "C"):
            gc_content += 1

    gc_content = gc_content / len(seq) * 100

    return gc_content


def is_acceptable_gc(seq: str, gc_bounds: tuple) -> bool:
    """
    Check if the GC content of a DNA sequence within specified bounds

    Parameters:
    - seq (str): DNA sequence
    - gc_bounds (tuple): tuple with lower and upper bounds for GC content

    Returns:
    - bool: True if GC content is within bounds, False otherwise
    """

    gc_content = calculate_gc_content(seq)
    lower_threshold, upper_threshold = gc_bounds

    res = lower_threshold <= gc_content <= upper_threshold

    return res


def is_acceptable_length(seq: str, length_bounds: Tuple[int, int]) -> bool:
    """
    Check if the length of a sequence falls within specified bounds.

    Parameters:
    - seq (str): A sequence.
    - length_bounds (tuple): A tuple representing the lower and upper bounds for sequence length.

    Returns:
    - bool: True if the length is within bounds, False otherwise.
    """

    read_len = len(seq)
    lower_threshold, upper_threshold = length_bounds

    res = lower_threshold <= read_len <= upper_threshold

    return res


encoding_dict = {
    "!": 0,
    '"': 1,
    "#": 2,
    "$": 3,
    "%": 4,
    "&": 5,
    "'": 6,
    "(": 7,
    ")": 8,
    "*": 9,
    "+": 10,
    ",": 11,
    "-": 12,
    ".": 13,
    "/": 14,
    "0": 15,
    "1": 16,
    "2": 17,
    "3": 18,
    "4": 19,
    "5": 20,
    "6": 21,
    "7": 22,
    "8": 23,
    "9": 24,
    ":": 25,
    ";": 26,
    "<": 27,
    "=": 28,
    ">": 29,
    "?": 30,
    "@": 31,
    "A": 32,
    "B": 33,
    "C": 34,
    "D": 35,
    "E": 36,
    "F": 37,
    "G": 38,
    "H": 39,
    "I": 40,
}


def is_acceptable_quality_score(seq: str, quality_threshold: int) -> bool:
    """
    Check if the average quality score of a sequence is above a specified threshold.

    Parameters:
    - seq (str): sequence of quality scores
    - quality_threshold (int): the minimum acceptable quality score

    Returns:
    - bool: True if the average quality score is above the threshold, False otherwise
    """

    read_quality_score = 0

    for score in seq:
        read_quality_score += encoding_dict[score]

    read_quality_score = read_quality_score / len(seq)

    res = read_quality_score >= quality_threshold

    return res


def write_fastq_file(
    fastq: Dict[str, Tuple[str, str]],
    input_path: str,
    output_filename: str = None,
):
    """
    Write a dictionary of sequences to a new FASTQ file.

    Parameters:
    - fastq (dict): Dictionary where keys are sequence identifiers and values are tuples containing sequence and quality scores.
    - input_path (str): Path to the input FASTQ file.
    - output_filename (str): Name of the output FASTQ file. If None, the base name of the input file is used.
    """

    if output_filename is None:
        output_filename = os.path.basename(input_path)

    data_dir = os.path.join(os.path.dirname(input_path), "fastq_filtrator_results")

    if not os.path.isdir(data_dir):
        os.mkdir(data_dir)

    with open(os.path.join(data_dir, output_filename), "w") as new_fastq:
        for key, item in fastq.items():
            new_fastq.write(key + "\n")
            new_fastq.write(item[0] + "\n")
            new_fastq.write("+" + key[1:] + "\n")
            new_fastq.write(item[1] + "\n")


def fastq_filter(
    input_path: str,
    output_filename: str = None,
    gc_bounds: tuple = (0, 100),
    length_bounds: tuple = (0, 2**32),
    quality_threshold: int = 0,
) -> dict:
    """
    Filter a dictionary of sequences based on specified criteria and write filtered dictionary in a FASTQ format.

    Parameters:
    - input_path (str): Path to the input FASTQ file
    - output_filename (str): Name of the output FASTQ file. If None, the base name of the input file is used
    - seqs (dict): dictionary where keys are sequence identifiers and values are tuples containing sequence and quality scores
    - gc_bounds (tuple): tuple with lower and upper bounds for GC content
    - length_bounds (tuple): tuple with lower and upper bounds for sequence length
    - quality_threshold (int): minimum acceptable quality score

    Returns:
    - dict: filtered dictionary containing only the sequences that meet the specified criteria
    """

    # fastq file to fastq dict
    fastq = read_fastq_file(input_path)

    # Check value of gc_bounds parameter
    if type(gc_bounds) == int:
        gc_bounds = (0, gc_bounds)

    elif len(gc_bounds) == 1:
        gc_bounds = (0, gc_bounds[0])

    elif len(gc_bounds) > 2:
        print("Error! Invalid gc_bound value!")
        return None

    # Check value of length_bounds parameter
    if type(length_bounds) == int:
        length_bounds = (0, length_bounds)

    elif len(length_bounds) == 1:
        length_bounds = (0, length_bounds[0])

    elif len(length_bounds) > 2:
        print("Error! Invalid length_bounds value!")
        return None

    # Check quality_threshold parameter
    if type(quality_threshold) != int:
        return "Error! Input integer value!"

    ids_to_del = []

    for id in fastq.keys():
        seq = fastq[id][0]
        quality = fastq[id][1]

        if (
            (not is_acceptable_gc(seq, gc_bounds))
            or (not is_acceptable_length(seq, length_bounds))
            or (not is_acceptable_quality_score(quality, quality_threshold))
        ):
            ids_to_del.append(id)
    for id in ids_to_del:
        del fastq[id]

    write_fastq_file(fastq, input_path, output_filename)

    return fastq
