import pandas as pd
import pytest
from deseq2_on_nf_RNAseq_functions import clean_gene_counts


# Define a fixture for a sample gene counts file
@pytest.fixture
def sample_gene_counts_file(tmp_path):
    # Create a temporary gene counts file
    file_path = tmp_path / "sample_gene_counts.tsv"
    with open(file_path, "w") as f:
        f.write(
            "gene_id\tgene_name\tsample1\tsample2\tsample3\n"
            "1\tGeneA\t0\t1\t0\n"
            "2\tGeneB\t0\t0\t0\n"
            "3\tGeneC\t5\t5\t5\n"
        )
    return str(file_path)


# Test the clean_gene_counts function
def test_clean_gene_counts(sample_gene_counts_file):
    # Call the function with test parameters
    gene_counts_df = clean_gene_counts(
        nextflow_gene_counts=sample_gene_counts_file,
        drop_unexpressed_genes=True,
        minimum_count_threshold=3,
    )

    # Verify the result by checking the DataFrame
    expected_result = pd.DataFrame({"GeneA": [1], "GeneC": [5]})

    # Compare the actual and expected DataFrames
    pd.testing.assert_frame_equal(gene_counts_df, expected_result)
