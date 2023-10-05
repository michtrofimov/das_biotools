def check_if_procedure(procedure):
    procedures = ["transcribe", "reverse", "complement", "reverse_complement"]

    if procedure not in procedures:
        raise ValueError(
            "This procedure is not supported yet! Check that you wrote your procedure correctly"
        )


def check_if_empty(seqs):
    # check if input is empty
    if len(seqs) == 0:
        raise ValueError("No sequences in input!")


def filter_seqs(seqs, dna_alphabet="ATGC", rna_alphabet="AUGC"):
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


def transcribe(seq):
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


def reverse_seq(seq):
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


def is_dna(seq):
    if ("T" in seq) or ("t" in seq):
        return True
    else:
        return False


def is_rna(seq):
    if ("U" in seq) or ("u" in seq):
        return True
    else:
        return False


def complement(seq):
    if (not is_dna(seq)) and (not is_rna(seq)):
        print("Ambigious sequence. Considering it as DNA")
        res = "".join([COMPLEMENTATION_DNA[nt] for nt in seq])

    elif is_dna(seq):
        res = "".join([COMPLEMENTATION_DNA[nt] for nt in seq])

    elif is_rna(seq):
        res = "".join([COMPLEMENTATION_RNA[nt] for nt in seq])

    return res


def reverse_complement(seq):
    res = complement(reverse_seq(seq))

    return res


def dna_rna_tools(*seqs):
    # parse seqs list
    *seqs, procedure = seqs

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
