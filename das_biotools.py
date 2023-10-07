from src.das_dna_rna_tools import dna_rna_tools
from src.das_protein_tools import protein_tools
from src.das_fastq_filter import fastq_filter

dna_rna_tools("atgc", "transcribe")
protein_tools("rqh", "get_pI")
fastq_filter(
    {
        "@SRX079804:1:SRR292678:1:1101:21885:21885": (
            "ACAGCAACATAAACATGATGGGATGGCGTAAGCCCCCGAGATATCAGTTTACCCAGGATAAGAGATTAAATTATGAGCAACATTATTAA",
            "FGGGFGGGFGGGFGDFGCEBB@CCDFDDFFFFBFFGFGEFDFFFF;D@DD>C@DDGGGDFGDGG?GFGFEGFGGEF@FDGGGFGFBGGD",
        ),
        "@SRX079804:1:SRR292678:1:1101:24563:24563": (
            "ATTAGCGAGGAGGAGTGCTGAGAAGATGTCGCCTACGCCGTTGAAATTCCCTTCAATCAGGGGGTACTGGAGGATACGAGTTTGTGTG",
            "BFFFFFFFB@B@A<@D>BDDACDDDEBEDEFFFBFFFEFFDFFF=CC@DDFD8FFFFFFF8/+.2,@7<<:?B/:<><-><@.A*C>D",
        ),
    },
    gc_bounds=40,
)
