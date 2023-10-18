from src.das_dna_rna_tools import dna_rna_tools
from src.das_protein_tools import protein_tools
from src.das_fastq_filter import fastq_filter

dna_rna_tools("atgc", procedure="complement")
protein_tools("rh", "rhqcqq", procedure="needleman_wunsch", gap_penalty=-2)
fastq_filter("example_fastq.fastq", gc_bounds=40)
