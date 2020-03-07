# Single-cell transcriptomics reveals new insights into the cell type-specific effects of obstructive (postvasectomy) and nonobstructive azoospermia in human testes

This repository contains scripts for data processing, analysis and figure generation using scRNA-Seq

## How to work with the data

These single-cell RNA sequencing data have been deposited on NODE project under the accession number [OEP000778](https://www.biosino.org/node/project/detail/OEP000778).

### 1. Obtaining the public scRNA-seq data
The paper of these datasets can be found at:

healthy: [Guo et al., Cell Rep. 2019](https://www.ncbi.nlm.nih.gov/pubmed/30315278)

obstructive azoospermia : [Sohni et al., Cell Rep. 2019](https://www.ncbi.nlm.nih.gov/pubmed/30726734)

##### Testes of Healthy individual (GSE112013)
```
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE112nnn/GSE112013/suppl/GSE112013_Combined_UMI_table.txt.gz 
```

##### Testes of OAa  (GSM3526587)
```
mkdir A1
cd A1
#GSM3526588_A1Total_barcodes.tsv.gz  
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3526nnn/GSM3526588/suppl/GSM3526588_A1Total_barcodes.tsv.gz
#GSM3526588_A1Total_genes.tsv.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3526nnn/GSM3526588/suppl/GSM3526588_A1Total_genes.tsv.gz
#GSM3526588_A1Total_matrix.mtx.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3526nnn/GSM3526588/suppl/GSM3526588_A1Total_matrix.mtx.gz
```
##### Testes of OAb  (GSM3526590)
```
mkdir A2
cd A2
#GSM3526588_A1Total_barcodes.tsv.gz  
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3526nnn/GSM3526590/suppl/GSM3526590_A2_total_barcodes.tsv.gz
#GSM3526588_A1Total_genes.tsv.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3526nnn/GSM3526590/suppl/GSM3526590_A2_total_genes.tsv.gz
#GSM3526588_A1Total_matrix.mtx.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM3526nnn/GSM3526590/suppl/GSM3526590_A2_total_matrix.mtx.gz
```
### 2. Data processing
To further process the data (normalization, dimensionality reduction and Integration), please follow either the [processing](../master/code/processing.r)
