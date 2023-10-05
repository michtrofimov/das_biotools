<img width="642" alt="Screenshot 2023-10-05 at 23 18 54" src="https://github.com/michtrofimov/das_biotools/assets/92677906/4410efa3-83fd-4a89-8ee7-6c6b6baa5805">

# DAS BIOTOOLS
> The great and terrifying successor of biopython

Das biotools is a homework project made during studies in Bioinformatics institute by Michil Trofimov in September 2023.

It is a python library which is able to perform several procedures on DNA/RNA, protein sequences and filter through FASTQ data.

# Features

## Tools for nucleotide sequences

- `check_if_procedure(procedure: str) -> None`: Checks if the specified procedure is supported. Raises a ValueError if the procedure is not supported.

- `check_if_empty(seqs: list) -> None`: Checks if the input sequence list is empty. Raises a ValueError if the input is empty.

- `filter_seqs(seqs: list, dna_alphabet: str = "ATGC", rna_alphabet: str = "AUGC") -> list`: Filters out sequences that are not nucleotide sequences. Returns a list of filtered nucleotide sequences.

- `transcribe(seq: str) -> str`: Transcribes a DNA sequence into an RNA sequence.

- `reverse_seq(seq: str) -> str`: Reverses a sequence.

- `complement(seq: str) -> str`: Generates the complement of a nucleotide sequence, considering it as either DNA or RNA.

- `reverse_complement(seq: str) -> str`: Generates the reverse complement of a nucleotide sequence, considering it as either DNA or RNA.

- `dna_rna_tools(*seqs: str, procedure: str) -> str`: Executes a specified procedure on a list of DNA/RNA sequences. The supported procedures are "transcribe," "reverse," "complement," and "reverse_complement."

### Usage
- **transcribe**
```python
dna_rna_tools('ATCG', procedure='transcribe') -> "AUCG"
```
- **reverse_seq**
```python
dna_rna_tools('ATCG', procedure='reverse') -> "GCTA"
```
- **complement**
```python
dna_rna_tools('ATCG', procedure='complement') -> "TAGC"
```

- **reverse_complement**
```python
dna_rna_tools('ATCG', procedure='reverse_complement') -> "CGAT"
```
## Tools for aminoacid sequences

- `is_protein(seq: str) -> bool`: Check if a given sequence is a valid protein sequence.
  
- `get_pI(sequence: str, pI_values: dict = None) -> str`: Gives isoelectric point value for each aminoacid individually.

- `build_scoring_matrix(match_score: int, mismatch_score: int, amino_acid_alphabet: str = None) -> dict`: Auxiliary function for needleman_wunsch. Build a scoring matrix for amino acid pairs, which can be used in sequence alignment algorithms.

- `needleman_wunsch(seq1: str, seq2: str) -> str`: Implement the Needleman-Wunsch algorithm for global sequence alignment of two amino acid sequences.

- `calculate_aa_freq(sequences: str) -> dict`: Calculate the frequences of aminoacids in protein sequences.

- `convert_to_3L_code(seq: str) -> str`: Converts one letter animoacid sequence to three letter aminoacid sequence.

- `protein_mass(seq: str) -> float`: Calculates molecular weight of the aminoacid sequence using monoisotopic masses.

- `translate_protein_rna(seq: str) -> str`: Converts aminoacid sequence to RNA sequence. For those aminoacids that are coded with more than one codon, this function randomly chooses one codon from the set.
- `protein_tools(*args: str) -> any`: Performs various actions on protein sequences. The supported procedures "get_pI", "needleman_wunsch", "build_scoring_matrix", "calculate_aa_freq", "translate_protein_rna", "convert_to_3L_code", "protein_mass".

### Usage

- **get_pI**

```python
protein_tools('RAAH', 'get_pI') -> "Sequence: RAAH. Isoelectric point of each aminoacid: [('R', 10.76), ('A', 6.0), ('H', 7.64)]"
```

- **build_scoring_matrix**
```python
protein_tools('2', '3', 'build_scoring_matrix') -> {'A': {'A': 2, 'C': 3, 'D': 3, ...}, 'C': {'A': 3, 'C': 2, 'D': 3, ...}, ...}
```

- **needleman_wunsch**
```python
protein_tools('RAAH', 'ARAH', 'needleman_wunsch') -> "RAAH, -ARA, final score: 1"
```

- **calculate_aa_freq**
```python
protein_tools('RAAH', 'calculate_aa_freq') -> {'R': 1, 'A': 2, 'H': 1}
```

- **convert_to_3L_code**
```python
protein_tools('RAAH', 'convert_to_3L_code') -> "ArgAlaAlaHis"
```

- **protein_mass**
```python
protein_tools('RAAH', 'protein_mass') -> 380.4934
```

- **translate_protein_rna**
```python
protein_tools('RAAH', 'translate_protein_rna') -> "GACGACGACAUGAC"
```

## Tools for FASTQ data filtering
- `calculate_gc_content(seq: str) -> float`: Calculates the GC content percentage of a DNA sequence.

- `is_acceptable_gc(seq: str, gc_bounds: tuple) -> bool`: Checks if the GC content of a DNA sequence is within specified bounds.

- `is_acceptable_length(seq: str, length_bounds: tuple) -> bool`: Checks if the length of a sequence falls within specified bounds.

- `is_acceptable_quality_score(seq: str, quality_threshold: int) -> bool`: Checks if the average quality score of a sequence is above a specified threshold.

- `fastq_filter(seqs: dict, gc_bounds: tuple = (0, 100), length_bounds: tuple = (0, 2**32), quality_threshold: int = 0) -> dict`: Filters a dictionary of sequences based on specified criteria. Returns a filtered dictionary containing only the sequences that meet the specified criteria.

### Usage 

- **calculate_gc_content**
```python
fastq_filter({'read1': ('ATCG', '!!@#'), 'read2': ('GCTA', '$$%&'), 'read3': ('AAAA', '%%%%')}, gc_bounds=(40, 60)) -> {'read1': ('ATCG', '!!@#'), 'read2': ('GCTA', '$$%&')}
```

- **is_acceptable_length**
```python
fastq_filter({'read1': ('ATCG', '!!@#'), 'read2': ('GCTA', '$$%&'), 'read3': ('AAAA', '%%%%')}, length_bounds=(2, 4)) -> {'read1': ('ATCG', '!!@#'), 'read3': ('AAAA', '%%%%')}
```

- **is_acceptable_quality_score**
```python
fastq_filter({'read1': ('ATCG', '!!@#'), 'read2': ('GCTA', '$$%&'), 'read3': ('AAAA', '%%%%')}, quality_threshold=10) -> {'read1': ('ATCG', '!!@#'), 'read2': ('GCTA', '$$%&'), 'read3': ('AAAA', '%%%%')}
```
