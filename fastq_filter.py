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


def is_acceptable_length(seq: str, length_bounds: tuple) -> bool:
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
