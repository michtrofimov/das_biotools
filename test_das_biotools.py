import pytest
import tempfile
import requests

from bio_files_processor import (
    convert_multiline_fasta_to_oneline,
    select_genes_from_gbk_to_fasta,
    OpenFasta,
    FastaRecord,
)
from das_biotools import DNASequence, RNASequence, AminoAcidSequence

@pytest.fixture
def input_data_fasta() -> list:
    """
    Fixture for input fasta data.
    """
    fasta = [
        FastaRecord(
            id=">GTD323452",
            description="5S_rRNA NODE_272_length_223_cov_0.720238:18-129(+)",
            seq="ACGGCCATAGGACTTTGAAAGCACCGCATCCCGTCCGATCTGCGAAGTTAACCAAGATGCCGCCTGGTTAGTACCATGGTGGGGGACCACATGGGAATCCCTGGTGCTGTG",
        ),
        FastaRecord(
            id=">GTD678345",
            description="16S_rRNA NODE_80_length_720_cov_1.094737:313-719(+)",
            seq="TTGGCTTCTTAGAGGGACTTTTGATGTTTAATCAAAGGAAGTTTGAGGCAATAACAGGTCTGTGATGCCCTTAGATGTTCTGGGCCGCACGCGCGCTACACTGAGCCCTTGGGAGTGGTCCATTTGAGCCGGCAACGGCACGTTTGGACTGCAAACTTGGGCAAACTTGGTCATTTAGAGGAAGTAAAAGTCGTAACAAGGT",
        ),
        FastaRecord(
            id=">GTD174893",
            description="16S_rRNA NODE_1_length_2558431_cov_75.185164:2153860-2155398(+)",
            seq="TTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGGGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT",
        ),
        FastaRecord(
            id=">GTD906783",
            description="16S_rRNA NODE_1_length_2558431_cov_75.185164:793941-795479(-)",
            seq="TTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGGGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT",
        ),
        FastaRecord(
            id=">GTD129563",
            description="16S_rRNA NODE_4_length_428221_cov_75.638017:281055-282593(-)",
            seq="CGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT",
        ),
    ]
    return fasta


def test_OpenFasta(input_data_fasta: list):
    """
    Test the OpenFasta class.
    """
    file_path = "./data/example_fasta.fasta"
    with OpenFasta(file_path) as fasta_file:
        records = fasta_file.read_records()

    assert len(records) == len(input_data_fasta)
    for record, expected_record in zip(records, input_data_fasta):
        assert record == expected_record


def test_run_genscan_correct_url():
    """
    Test if ValueError raises when the Genscan URL is not accessible.
    """
    url = "http://argonaute.mit.edu/cgi-bin/genscanw_py.cgi_WRONG_URL"
    with pytest.raises(ValueError):
        response = requests.post(url)
        status = response.status_code

        if status != 200:
            raise ValueError("Failed to connect to GENSCAN service")
    # assert response.status_code == 200


@pytest.fixture
def gbk_file(tmpdir):
    """
    Fixture for creating a temporary GBK file.
    """
    gbk_path = tmpdir.join("example_gbk.gbk")
    with open(gbk_path, "w") as gbk:
        gbk.write(
            ">Feature sample\n"
            "    CDS             123..456\n"
            '                    /gene="first"\n'
            '                    /translation="ATGCAT"\n'        
            ">Feature sample\n"
            "    CDS             123..456\n"
            '                    /gene="sample"\n'
            '                    /translation="ATGC"\n'
            ">Feature sample2\n"
            "    CDS             789..1012\n"
            '                    /gene="test"\n'
            '                    /translation="ATGCCCC"\n'
        )
    return gbk_path


def test_select_genes_from_gbk_to_fasta_output_len(gbk_file):
    """
    Test the output length of select_genes_from_gbk_to_fasta function.
    """
    output_file = tempfile.NamedTemporaryFile(delete=False).name
    select_genes_from_gbk_to_fasta(
        str(gbk_file), ["sample"], output_fasta=output_file
    )  # Convert Path object to string
    with open(output_file, "r") as output_fasta:
        lines = output_fasta.readlines()
        assert len(lines) == 4

def test_select_genes_from_gbk_to_fasta_flanking_genes_names(gbk_file):
    """
    Test the flanking genes names in the output of select_genes_from_gbk_to_fasta function.
    """
    output_file = tempfile.NamedTemporaryFile(delete=False).name
    select_genes_from_gbk_to_fasta(
        str(gbk_file), "sample", output_fasta=output_file
    )  # Convert Path object to string
    with open(output_file, "r") as output_fasta:
        lines = output_fasta.readlines()
        assert lines[0] == ">first\n"
        assert lines[2] == ">test\n"

def test_select_genes_from_gbk_to_fasta_output_gene_sequence(gbk_file):
    """
    Test the output gene sequence of select_genes_from_gbk_to_fasta function.
    """
    output_file = tempfile.NamedTemporaryFile(delete=False).name
    select_genes_from_gbk_to_fasta(
        str(gbk_file), "sample", output_fasta=output_file
    )  # Convert Path object to string
    with open(output_file, "r") as output_fasta:
        lines = output_fasta.readlines()
        assert lines[1] == "ATGCAT\n"
        assert lines[3] == "ATGCCCC\n"
        


def test_transcribe():
    """
    Test transcribe function of DNASequence class.
    """
    sequence = "ATGC"
    rna_seq = "AUGC"
    dna_seq = DNASequence(sequence)
    expected_transcribed_seq = RNASequence(rna_seq)
    transcribed_seq = dna_seq.transcribe()
    assert expected_transcribed_seq.sequence == transcribed_seq.sequence


def test_alphabet_checking():
    """
    Test alphabet_checking method of RNASequence class.
    """
    sequence = "AUGCU"
    rna_seq = RNASequence(sequence)
    assert rna_seq.is_valid_alphabet() == True

def test_calculate_aa_freq():
    """
    Test alphabet_checking method of RNASequence class.
    """
    sequence = "LIMMM"
    aa_seq = AminoAcidSequence(sequence)
    aa_freq = aa_seq.calculate_aa_freq()
    aa_freq_M = aa_freq['M']
    assert aa_freq_M == 3
