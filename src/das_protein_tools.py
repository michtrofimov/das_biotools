# importing necessary modules
import src.protein_dict as protein_dict
from random import choice

AMINO_LETTERS = set("ACDEFGHIKLMNPQRSTVWY")


# Function to determine is the sequence is a protein or not
def is_protein(seq: str) -> bool:
    """
    This function checks if the sequence is a protein or not

    Arguments:
        seq (str): A sequence of aminoacids

    Output:
        returns True or False
    """
    unique_chars = set(seq)
    return unique_chars.issubset(AMINO_LETTERS)


# Function to get pI for each aa
def get_pI(
    sequence: str,
    pI_values: dict = None,
) -> str:
    """
    Gives isoelectric point value for each aminoacid individually

    Args:
    - sequence (str): sequence for which to calculate isoelectric point
    - pI_values (dict): acid dissociation constants for each aminoacid

    Return:
    - str: string, which contains:
            - an original sequence,
            - list of tuple pairs of aminoacid and corresponding isoelectric point,
    """

    if pI_values is None:
        # Default pKa_values if not provided
        pI_values = protein_dict.aa_pI

    aminoacid_pIs = []

    # Calculate pI for each amino acid in the sequence while preserving case
    analysed_aa = []
    for aa in sequence:
        aa_upper = aa.upper()
        if aa_upper not in analysed_aa:
            if aa_upper in pI_values:
                pI = pI_values[aa_upper]
                analysed_aa.append(aa_upper)
                if aa.isupper():
                    aminoacid_pIs.append((aa_upper, pI))
                else:
                    aminoacid_pIs.append((aa, pI))
        else:
            continue

    return f"Sequence: {sequence}. Isoelectric point of each aminoacid: {aminoacid_pIs}"


# Function to build scoring matrix for needleman_wunsch function
def build_scoring_matrix(
    match_score: int,
    mismatch_score: int,
    amino_acid_alphabet: str = None,
) -> dict:
    """
    Build a default scoring matrix, if not provided in needleman-wunsch function parameter

    Args:
    - match_score (int): integer value of a matching score of aminoacids
    - mismatch_score (int): integer value of a mismatching score of aminoacids
    - amino_acid_alphabet (str): upper case amino acid alphabet

    Returns:
    - a dictionary of dictionaries representing a scoring matrix for aminoacids paris. Key of a dictionary is an aminoacid and its value is a dictionary of scores
    """

    if amino_acid_alphabet is None:
        amino_acid_alphabet = AMINO_LETTERS

    scoring_matrix = {}

    for aa1 in amino_acid_alphabet:
        scoring_matrix[aa1] = {}
        for aa2 in amino_acid_alphabet:
            scoring_matrix[aa1][aa2] = (
                match_score if aa1.upper() == aa2.upper() else mismatch_score
            )

    return scoring_matrix


# Function to perform alignment based on needleman_wunsch algorithm
def needleman_wunsch(
    seq1: str,
    seq2: str,
    scoring_matrix: dict = None,
    gap_penalty: int = -1,
    match_score: int = 1,
    mismatch_score: int = -1,
) -> str:
    """
    Uses Needleman-Wunsch algorithm to make a global alignment of two sequences.

    Args:
    - seq1 (str): first aminoacid sequence for alignment
    - seq2 (str): second aminoacid sequence for alignment
    - scoring_matrix (dict): A dictionary representing a scoring matrix for amino acid pairs
      If not provided, a default scoring_matrix is generated based on match and mismatch scores
    - gap_penalty (int): integer va;ue of a penalty score for introducing a gap in the alignment
    - match_score (int): integer value of a matching score for matching aminoacids
    - mismatch_score (int): integer value of a mismatching score for mismatched aminoacids

    Returns:
    - string: a string containing the aligned sequences (str), the aligned score (int)
    """
    if scoring_matrix is None:
        # Default scoring matrix if not provided
        scoring_matrix = build_scoring_matrix(match_score, mismatch_score)

    seq1_upper = seq1.upper()  # Convert seq1 to uppercase
    seq2_upper = seq2.upper()  # Convert seq2 to uppercase

    m, n = len(seq1_upper), len(seq2_upper)

    # Initialize matrices
    dp = [[0] * (n + 1) for _ in range(m + 1)]
    traceback = [[""] * (n + 1) for _ in range(m + 1)]

    # Fill in the scoring matrix and traceback matrix
    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = dp[i - 1][j - 1] + scoring_matrix.get(seq1_upper[i - 1], {}).get(
                seq2_upper[j - 1], mismatch_score
            )
            delete = dp[i - 1][j] + gap_penalty
            insert = dp[i][j - 1] + gap_penalty

            dp[i][j] = max(match, delete, insert)

            if dp[i][j] == match:
                traceback[i][j] = "D"  # Diagonal (indicates a match/mismatch)
            elif dp[i][j] == delete:
                traceback[i][j] = "U"  # Up (indicates a gap in seq2)
            else:
                traceback[i][j] = "L"  # Left (indicates a gap in seq1)

    # Traceback to find the aligned sequences while preserving case
    aligned_seq1, aligned_seq2 = [], []
    i, j = m, n
    while i > 0 or j > 0:
        if traceback[i][j] == "D":
            # check original case of amionacid in seq1
            if seq1[i - 1].isupper():
                aligned_seq1.append(seq1_upper[i - 1])
            else:
                aligned_seq1.append(seq1[i - 1])

            # check original case of amionacid in seq2
            if seq2[j - 1].isupper():
                aligned_seq2.append(seq2_upper[j - 1])
            else:
                aligned_seq2.append(seq2[j - 1])

            i -= 1
            j -= 1
        elif traceback[i][j] == "U":
            # check original case of amionacid in seq1
            if seq1[i - 1].isupper():
                aligned_seq1.append(seq1_upper[i - 1])
            else:
                aligned_seq1.append(seq1[i - 1])
            aligned_seq2.append("-")

            i -= 1
        else:
            aligned_seq1.append("-")
            # check original case of amionacid in seq2
            if seq2[j - 1].isupper():
                aligned_seq2.append(seq2_upper[j - 1])
            else:
                aligned_seq2.append(seq2[j - 1])

            j -= 1

    aligned_seq1 = "".join(reversed(aligned_seq1))
    aligned_seq2 = "".join(reversed(aligned_seq2))

    return f"{aligned_seq1}, {aligned_seq2}, final score: {dp[m][n]}"


# Function to calculate frequency of unique aminoacid in the sequence
def calculate_aa_freq(sequences: str) -> dict:
    """
    Calculates the frequency of each amino acid in a protein sequence or sequences.

    :param sequences: protein sequence or sequences
    :type sequences: str or list of str
    :return: dictionary with the frequency of each amino acid
    :rtype: dict
    """

    # Creating a dictionary with aminoacid frequencies:
    amino_acid_frequency = {}

    for amino_acid in sequences:
        # If the aminoacid has been already in:
        if amino_acid in amino_acid_frequency:
            amino_acid_frequency[amino_acid] += 1
        # If the aminoacid hasn't been already in:
        else:
            amino_acid_frequency[amino_acid] = 1

    return amino_acid_frequency


# Function to convert one-letter protein sequence to three-letter protein sequence
def convert_to_3L_code(seq: str) -> str:
    """
    This function takes one letter aminoacids sequence and convert's it to three leter coding

    Arguments:
        seq (str): A sequence of aminoacids

    Output:
        same sequence but in three-letter coding
    """
    seq = seq.upper()
    if is_protein(seq) is True:
        sequence = "".join(protein_dict.aa_one_to_three_letter.get(aa) for aa in seq)
        return sequence[:-1]
    else:
        raise ValueError("Sequence is not a protein, input should be protein")


# Function to calculate protein mass
def protein_mass(seq: str) -> float:
    """
    This function takes aminoacids sequence and counts it's summary molecular weight using monoisotopic masses

    Arguments:
        seq (str): A sequence of aminoacids

    Output:
        returns molecular weight
    """
    seq = seq.upper()
    if is_protein(seq) is True:
        mass = sum(protein_dict.aa_monoistopic_mass_dict.get(aa) for aa in seq)
        return mass
    else:
        raise ValueError("Sequence is not a protein, input should be protein")


# Function to translate Protein to RNA
def translate_protein_rna(seq: str) -> str:
    """
    This function takes  aminoacid sequence and translates in to the RNA.
    As most of the aminoacids are coded with several different codons,
    this function will take a random codon of the set for such aminoacids.

    Arguments:
        seq (str): A sequence of RNA molecule

    Output:
        returns sequence of aminoacids
    """
    seq = seq.upper()
    if is_protein(seq) is True:
        rna = ""
        for aa in seq:
            codon = choice(protein_dict.aa_codon_dict.get(aa))
            rna += codon
        return rna
    else:
        raise ValueError("Sequence is not a protein, input should be a protein")


def protein_tools(*args: any, procedure: str):
    """
    Main function to perform various procedures on protein sequences.

    Args:
    - *args: Variable number of arguments. The first one or two arguments should be protein sequences,
             other arguments are auxiliary arguments.
    - procedure (str):  A procedure to be performed
    Returns:
    - The result of the specified procedure on the input protein sequences.

    Raises:
    - ValueError: If the specified procedure is not supported or if there is an error in the number of sequences.
                  Also raised if the input sequences are not valid protein sequences.

    Supported procedures:
    - "get_pI": Calculate isoelectric points for each amino acid in the sequence.
    - "needleman_wunsch": Perform global alignment of two sequences using the Needleman-Wunsch algorithm.
    - "build_scoring_matrix": Build a scoring matrix for amino acid pairs.
    - "calculate_aa_freq": Calculate the frequency of each amino acid in a protein sequence.
    - "translate_protein_rna": Translate amino acid sequence to RNA, using random codons for each amino acid.
    - "convert_to_3L_code": Convert one-letter amino acid sequence to three-letter coding.
    - "protein_mass": Calculate the molecular weight of the protein sequence.
    """

    procedure_list = {
        "get_pI": get_pI,
        "needleman_wunsch": needleman_wunsch,
        "build_scoring_matrix": build_scoring_matrix,
        "calculate_aa_freq": calculate_aa_freq,
        "translate_protein_rna": translate_protein_rna,
        "convert_to_3L_code": convert_to_3L_code,
        "protein_mass": protein_mass,
    }

    # Check whether procedure is valid
    if procedure not in procedure_list:
        raise ValueError(f"No such action: {procedure}")

    # Check whether number of args for functions is valid
    if not (procedure == "needleman_wunsch" and len(args) < 2):
        raise ValueError("Error in number of sequences")
    elif not (procedure == "get_pI" and len(args) > 3):
        raise ValueError("Error in number of sequences")
    elif not (
        procedure != "needleman_wunsch" and procedure != "get_pI" and len(args) == 2
    ):
        raise ValueError("Error in number of sequences")

    # Check whether query sequence is protein
    if not (procedure == "needleman_wunsch"):
        if not is_protein(args[0]):
            raise ValueError(f"The sequence is not protein sequence: {args[0]}")

        result = procedure_list[procedure](*args)
        return result

    else:
        if not (is_protein(*args[0]) or is_protein(*args[1])):
            raise ValueError(f"Sequences are not protein sequences")

        result = procedure_list[procedure](*args)

    return result
