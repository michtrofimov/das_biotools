def check_if_procedure(procedure: str) -> None:
    """
    Check if the specified procedure is valid

    Parameters:
    - procedure (str): DNA/RNA processing procedure

    Raises:
    - ValueError: if the procedure is not in the supported list
    """
    procedures = ["transcribe", "reverse", "complement", "reverse_complement"]

    if procedure not in procedures:
        raise ValueError(
            "This procedure is not supported yet! Check that you wrote your procedure correctly"
        )


def check_if_empty(seqs: list) -> None:
    """
    Check if the input sequence list is empty

    Parameters:
    - seqs (list): DNA/RNA sequences

    Raises:
    - ValueError: if the input sequence list is empty
    """

    # check if input is empty
    if len(seqs) == 0:
        raise ValueError("No sequences in input!")


def filter_seqs(
    seqs: list, dna_alphabet: str = "ATGC", rna_alphabet: str = "AUGC"
) -> list:
    """
    Filter out non-nucleotide sequences from the input list

    Parameters:
    - seqs (list): DNA/RNA sequences
    - dna_alphabet (str): DNA alphabet characters
    - rna_alphabet (str): RNA alphabet characters

    Returns:
    - list: filtered list of DNA/RNA sequences
    """

    check_if_empty(seqs)
    for seq in seqs:
        # check if sequence is not a nucleotide sequence and pop it
        seq_letters = set(seq.upper())
        is_dna = seq_letters.issubset(dna_alphabet)
        is_rna = seq_letters.issubset(rna_alphabet)

        if (not is_dna) and (not is_rna):
            print(f"Sequence {seq} is not a nucleotide sequence!")
            seqs.remove(seq)
    check_if_empty(seqs)

    return seqs


def transcribe(seq: str) -> str:
    """
    transcribe a DNA sequence to an RNA sequence

    Parameters:
    - seq (str): DNA sequence

    Returns:
    - str: transcribed RNA sequence
    """

    transcribed = ""

    for nt in seq:
        if nt == "T":
            transcribed += "U"

        elif nt == "t":
            transcribed += "u"

        elif (nt == "U") or (nt == "u"):
            return "Error! RNA sequence for transcription."

        else:
            transcribed += nt

    return transcribed


def reverse_seq(seq: str) -> str:
    """
    Reverse a DNA/RNA sequence

    Parameters:
    - seq (str): input DNA/RNA sequence

    Returns:
    - str: reversed sequence
    """

    return seq[::-1]


COMPLEMENTATION_DNA = {
    "A": "T",
    "G": "C",
    "T": "A",
    "C": "G",
    "a": "t",
    "g": "c",
    "t": "a",
    "c": "g",
}
COMPLEMENTATION_RNA = {
    "A": "U",
    "G": "C",
    "U": "A",
    "C": "G",
    "a": "u",
    "g": "c",
    "u": "a",
    "c": "g",
}


def is_dna(seq: str) -> bool:
    """
    Check if a sequence is DNA

    Parameters:
    - seq (str): input sequence

    Returns:
    - bool: True if the sequence is DNA, False otherwise
    """

    if ("T" in seq) or ("t" in seq):
        return True
    else:
        return False


def is_rna(seq: str) -> bool:
    """
    Check if a sequence is RNA

    Parameters:
    - seq (str): Input sequence

    Returns:
    - bool: True if the sequence is RNA, False otherwise
    """

    if ("U" in seq) or ("u" in seq):
        return True
    else:
        return False


def complement(seq: str) -> str:
    """
    Generate the complement of a DNA/RNA sequence

    Parameters:
    - seq (str): Input sequence.

    Returns:
    - str: Complemented sequence.
    """

    if (not is_dna(seq)) and (not is_rna(seq)):
        print("Ambigious sequence. Considering it as DNA")
        res = "".join([COMPLEMENTATION_DNA[nt] for nt in seq])

    elif is_dna(seq):
        res = "".join([COMPLEMENTATION_DNA[nt] for nt in seq])

    elif is_rna(seq):
        res = "".join([COMPLEMENTATION_RNA[nt] for nt in seq])

    return res


def reverse_complement(seq: str) -> str:
    """
    Generate the reverse complement of a DNA/RNA sequence

    Parameters:
    - seq (str): Input sequence

    Returns:
    - str: Reverse complemented sequence
    """

    res = complement(reverse_seq(seq))

    return res


def dna_rna_tools(*seqs: str, procedure: str) -> str:
    """
    Perform DNA/RNA operations based on the specified procedure

    Parameters:
    - *seqs (str): Variable number of DNA/RNA sequences
    - procedure (str): A DNA/RNA processing procedure

    Returns:
    - str: Result of the specified procedure on the sequences.

    Supported Procedures:
    - "transcribe": Transcribe DNA sequences to RNA
    - "reverse": Reverse the sequences
    - "complement": Generate the complement of DNA/RNA sequences
    - "reverse_complement": Generate the reverse complement of DNA/RNA sequences

    Raises:
    - ValueError: If the specified procedure is not supported.
    """

    # check whether a valid procedure
    check_if_procedure(procedure)

    # filter out not nucleotide sequences
    seqs = filter_seqs(seqs)

    procedures = {
        "transcribe": transcribe,
        "reverse": reverse_seq,
        "complement": complement,
        "reverse_complement": reverse_complement,
    }

    for seq in seqs:
        res = procedures[procedure](seq)

    return res
