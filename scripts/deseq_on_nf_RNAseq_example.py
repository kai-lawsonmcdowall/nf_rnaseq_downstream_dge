# %%
import pandas as pd
from deseq2_on_nf_RNAseq_functions import clean_gene_counts, deseq2_dataframe_generator

# %%
ipsc_macrophage_gene_counts = clean_gene_counts(
    "nf_RNAseq_outputs/salmon.merged.gene_counts_length_scaled.tsv",
)
ipsc_macrophage_gene_counts

# %%
ipsc_macrophage_experimental_design = pd.read_csv(
    "nf_RNAseq_outputs/test_clinical.tsv", sep="\t", header=0, index_col="SRR"
)
ipsc_macrophage_experimental_design

# %%
result = deseq2_dataframe_generator(
    cleaned_gene_counts_df=ipsc_macrophage_gene_counts,
    metadata=ipsc_macrophage_experimental_design,
    metadata_grouping="Description",
)
result

# %%
