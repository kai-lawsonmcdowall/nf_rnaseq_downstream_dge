import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats


def clean_gene_counts(
    nextflow_gene_counts: str,
    drop_unexpressed_genes: bool = True,
    minimum_count_threshold: int = None,
):
    """
    input:
        salmon_merged.gene_counts.tsv from Nextflow Pipeline
    output:
        gene counts df cleaned for downstream use in PyDESeq2.
    """

    # Read the gene counts file into a DataFrame
    gene_counts = pd.read_csv(nextflow_gene_counts, sep="\t", header=0, index_col=False)

    # Convert numeric columns to integers (except "gene_id" and "gene_name")
    numeric_cols = gene_counts.columns.difference(["gene_id", "gene_name"])
    gene_counts[numeric_cols] = gene_counts[numeric_cols].astype(int)

    # Drop unexpressed genes if specified
    if drop_unexpressed_genes:
        gene_counts = gene_counts[gene_counts[numeric_cols].any(axis=1)]

    # Apply minimum count threshold if specified
    if minimum_count_threshold:
        gene_counts = gene_counts[
            gene_counts[numeric_cols].sum(axis=1) >= minimum_count_threshold
        ]

    # Pivot the DataFrame for gene_names and gene_id columns
    gene_counts = gene_counts.pivot_table(columns="gene_name").astype(int)

    return gene_counts


def deseq2_dataframe_generator(
    cleaned_gene_counts_df: pd.DataFrame,
    metadata: pd.DataFrame,
    metadata_grouping: str,
):
    """f(x) to in-take cleaned counts and metadata and apply pydeseq2 to it.

    input:
        cleaned_gene_counts_df (pd.DataFrame): The counts matrix processed in clean_gene_nextflow_counts.
        clinical (pd.DataFrame): the metadata
        metadata_grouping (str): which column to make the comparison on. For example, you could have column x, which has control and samples.
    output
        stats_res (DeseqStatsobject )
    """
    # create DEseqDataSet object
    dds = DeseqDataSet(
        counts=cleaned_gene_counts_df,
        metadata=metadata,
        design_factors=metadata_grouping,  # the comparator.
    )

    # run deseq2 method to fit dispersions and LFCs.
    dds.deseq2()

    # compute p-values and adjusted p-values for differential expresion.
    stat_res = DeseqStats(dds)
    summary_df = stat_res.summary()

    # returns the
    return summary_df
