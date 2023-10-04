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
