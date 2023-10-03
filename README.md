### RNAseq pipeline outputs to DGE by pyDeSeq2 

The outputs of the nextflow RNAseq pipeline which we expect for downstream analysis. In particular here: 

- pyDeseq2


We will have to manually construct the metadata csv for input into pydeseq2 however, as the pipeline only provides the raw counts. 
`The deseq2_on_nf_RNAseq_outputs.py` expects a raw counts file, such as: 

- `/pathway_analysis/nf_RNAseq_outputs/salmon.merged.gene_counts_length_scaled.tsv` which is taken from the IPSC macrophage repo.

We could potentially use the salmon_tx2gene for additional downstream annotation. 


### To-Do:
- create a small test file from the nf_RNAseq_inputs/salmon.merged.gene_counts_length_scaled.tsv file 
- feed this into the test_gene_counts.py test
- make sure this test is accessing this file rather than generating a new one. 

