<img width="642" alt="Screenshot 2023-10-05 at 23 18 54" src="https://github.com/michtrofimov/das_biotools/assets/92677906/4410efa3-83fd-4a89-8ee7-6c6b6baa5805">

# DAS BIOTOOLS
"The great and terrifying successor of biopython"

Das biotools is a homework project made during studies in Bioinformatics institute by Michil Trofimov in September 2023.

It is a python library which is able to perform several procedures on DNA/RNA, protein sequences, filter through FASTQ data, convert multi-line FASTA to single-line FASTA and extract features from .gbk file.

# Installation
```bash
git clone git@github.com:michtrofimov/das_biotools.git
```

```bash
cd das_biotools
```
# Features

Tools for nucleotide sequences and aminoacids sequences support upper and lower cases.

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
  
- `get_pI(sequence: str, pI_values: dict = None) -> str`: Gives isoelectric point value for each aminoacid individually. User can pass their own pI values

- `build_scoring_matrix(match_score: int, mismatch_score: int, amino_acid_alphabet: str = None) -> dict`: Auxiliary function for needleman_wunsch. Build a scoring matrix for amino acid pairs, which can be used in sequence alignment algorithms.

- `needleman_wunsch(seq1: str,
    seq2: str,
    scoring_matrix: dict = None,
    gap_penalty: int = -1,
    match_score: int = 1,
    mismatch_score: int = -1) -> str`: Implement the Needleman-Wunsch algorithm for global sequence alignment of two amino acid sequences.

- `calculate_aa_freq(sequences: str) -> dict`: Calculate the frequences of aminoacids in protein sequences.

- `convert_to_3L_code(seq: str) -> str`: Converts one letter animoacid sequence to three letter aminoacid sequence.

- `protein_mass(seq: str) -> float`: Calculates molecular weight of the aminoacid sequence using monoisotopic masses.

- `translate_protein_rna(seq: str) -> str`: Converts aminoacid sequence to RNA sequence. For those aminoacids that are coded with more than one codon, this function randomly chooses one codon from the set.
- `protein_tools(*args: any, **kwargs: any) -> any`: Performs various actions on protein sequences.
    - *args: Variable number of arguments. The first one or two arguments should be protein sequences
    - **kwargs:  Keyword argumets. First one is procedure, other are additional arguments for specific functions
    - The supported procedures "get_pI", "needleman_wunsch", "build_scoring_matrix", "calculate_aa_freq", "translate_protein_rna", "convert_to_3L_code", "protein_mass"

### Usage

- **get_pI**

```python
protein_tools('rh', pI_values = {'R' : 4,'H' : 3.5}, procedure='get_pI') -> "Sequence: rh. Isoelectric point of each aminoacid: [('r', 4), ('h', 3.5)]:
```

- **build_scoring_matrix**
```python
protein_tools('2', '3', procedure='build_scoring_matrix') -> {'A': {'A': 2, 'C': 3, 'D': 3, ...}, 'C': {'A': 3, 'C': 2, 'D': 3, ...}, ...}
```

- **needleman_wunsch**
```python
protein_tools('rh','rhqcqq',procedure='needleman_wunsch', gap_penalty = -2, match_score = 2) -> '----rh, rhqcqq, final score: -2'
```

- **calculate_aa_freq**
```python
protein_tools('RAAH', procedure='calculate_aa_freq') -> {'R': 1, 'A': 2, 'H': 1}
```

- **convert_to_3L_code**
```python
protein_tools('RAAH', procedure='convert_to_3L_code') -> "ArgAlaAlaHis"
```

- **protein_mass**
```python
protein_tools('RAAH', procedure='protein_mass') -> 380.4934
```

- **translate_protein_rna**
```python
protein_tools('RAAH', procedure='translate_protein_rna') -> "GACGACGACAUGAC"
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

## Tools for bio-files processing

- `convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None)`: Convert a multi-line FASTA file to a one-line FASTA file.

- `def select_genes_from_gbk_to_fasta(
    input_gbk: str,
    genes: List[str],
    n_before: int = 1,
    n_after: int = 1,
    output_fasta: str = None,
):`: Extracts gene sequences from a GenBank (gbk) file and creates a FASTA file with specified neighboring genes.

### Usage

- **convert_multiline_fasta_to_oneline**
```python
convert_multiline_fasta_to_oneline('example_multiline_fasta.fasta') -> "example_multiline_fasta_converted.fasta"
```

- **select_genes_from_gbk_to_fasta**
```python
select_genes_from_gbk_to_fasta('example_gbk.gbk', 'flu_2', n_before = 2, n_after = 3) -> "example_gbk_flanking_genes.fasta"
```
