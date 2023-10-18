# Dictionary: keys - single-letter amino acid designations; values - list of RNA codons
aa_codon_dict = {
    "G": ["GGA", "GGU", "GGC", "GGG"],
    "R": ["AGA", "AGG", "CGA", "CGC", "CGG", "CGU"],
    "S": ["AGC", "AGU", "UCA", "UCC", "UCG", "UCU"],
    "E": ["GAA", "GAG"],
    "P": ["CCA", "CCC", "CCG", "CCU"],
    "L": ["CUA", "CUC", "CUG", "CUU", "UUA", "UUG"],
    "V": ["GUA", "GUC", "GUG", "GUU"],
    "T": ["ACA", "ACC", "ACG", "ACU"],
    "A": ["GCA", "GCC", "GCG", "GCU"],
    "I": ["AUA", "AUC", "AUU"],
    "F": ["UUC", "UUU"],
    "H": ["CAC", "CAU"],
    "Y": ["UAC", "UAU"],
    "Q": ["CAA", "CAG"],
    "C": ["UGC", "UGU"],
    "N": ["AAC", "AAU"],
    "D": ["GAC", "GAU"],
    "K": ["AAA", "AAG"],
    "M": ["AUG"],
    "W": ["UGG"],
}


# Dictionary: keys - single-letter amino acid designations; values - names of amino acids
aa_one_to_three_letter = {
    "A": "Ala-",
    "C": "Cys-",
    "D": "Asp-",
    "E": "Glu-",
    "F": "Phe-",
    "G": "Gly-",
    "H": "His-",
    "I": "Ile-",
    "K": "Lys-",
    "L": "Leu-",
    "M": "Met-",
    "N": "Asn-",
    "P": "Pro-",
    "Q": "Gln-",
    "R": "Arg-",
    "S": "Ser-",
    "T": "Thr-",
    "V": "Val-",
    "W": "Trp-",
    "Y": "Tyr-",
}


# aminoacids mass dictionary
aa_monoistopic_mass_dict = {
    "A": 71.03711,
    "C": 103.00919,
    "D": 115.02694,
    "E": 129.04259,
    "F": 147.06841,
    "G": 57.02146,
    "H": 137.05891,
    "I": 113.08406,
    "K": 128.09496,
    "L": 113.08406,
    "M": 131.04049,
    "N": 114.04293,
    "P": 97.05276,
    "Q": 128.05858,
    "R": 156.10111,
    "S": 87.03203,
    "T": 101.04768,
    "V": 99.06841,
    "W": 186.07931,
    "Y": 163.06333,
}

# aminoacids pI (isoelectric point) values dictionary
aa_pI = {
    "A": 6.0,  # Alanine
    "R": 10.8,  # Arginine
    "N": 5.4,  # Asparagine
    "D": 2.8,  # Aspartic Acid
    "C": 5.0,  # Cysteine
    "E": 3.2,  # Glutamic Acid
    "Q": 5.7,  # Glutamine
    "G": 6.1,  # Glycine
    "H": 7.6,  # Histidine
    "I": 6.0,  # Isoleucine
    "L": 6.0,  # Leucine
    "K": 9.7,  # Lysine
    "M": 5.7,  # Methionine
    "F": 5.5,  # Phenylalanine
    "P": 6.3,  # Proline
    "S": 5.7,  # Serine
    "T": 5.6,  # Threonine
    "W": 5.9,  # Tryptophan
    "Y": 5.7,  # Tyrosine
    "V": 6.0,  # Valine
}
