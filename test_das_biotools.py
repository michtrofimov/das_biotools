import pytest
import tempfile
import requests

from pathlib import Path
from bio_files_processor import (
    convert_multiline_fasta_to_oneline,
    select_genes_from_gbk_to_fasta,
    OpenFasta,
    FastaRecord,
)

from das_biotools import DNASequence, RNASequence, AminoAcidSequence


@pytest.fixture
def input_data_fasta():
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


def test_OpenFasta(input_data_fasta):
    file_path = "./data/example_fasta.fasta"
    with OpenFasta(file_path) as fasta_file:
        records = fasta_file.read_records()

    assert len(records) == len(input_data_fasta)
    for record, expected_record in zip(records, input_data_fasta):
        assert record == expected_record


def test_run_genscan_correct_url():
    """
    Test if the Genscan URL is accessible.
    """
    url = "http://argonaute.mit.edu/cgi-bin/genscanw_py.cgi"
    response = requests.post(url)

    assert response.status_code == 200


@pytest.fixture
def gbk_file(tmp_path):
    gbk_path = tmp_path / "test.gbk"
    with open(gbk_path, "w") as gbk:
        gbk.write(
            ">Feature sample\n"
            "    CDS             123..456\n"
            '                    /gene="sample"\n'
            '                    /translation="ATGC"\n'
            ">Feature sample2\n"
            "    CDS             789..1012\n"
            '                    /gene="sample2"\n'
            '                    /translation="ATGC"\n'
        )
    return gbk_path


def test_select_genes_from_gbk_to_fasta(gbk_file):
    output_file = tempfile.NamedTemporaryFile(delete=False).name
    select_genes_from_gbk_to_fasta(
        str(gbk_file), ["sample"], output_fasta=output_file
    )  # Convert Path object to string
    with open(output_file, "r") as output_fasta:
        lines = output_fasta.readlines()
        assert len(lines) == 2
        assert lines[0] == ">sample\nATGC\n"
        assert lines[1] == ">sample2\nATGC\n"


def test_transcribe():
    """
    Test transcribe function of DNASequence class.
    """
    sequence = "ATGC"
    dna_seq = DNASequence(sequence)
    expected_transcribed_seq = RNASequence("AUGC")
    transcribed_seq = dna_seq.transcribe()
    assert expected_transcribed_seq == transcribed_seq


def test_alphabet_checking():
    """
    Test alphabet_checking method of RNASequence class.
    """
    sequence = "AUGCU"
    rna_seq = RNASequence(sequence)
    assert rna_seq.is_valid_alphabet() == True
