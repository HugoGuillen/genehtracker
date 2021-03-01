# genehtrack - Gene ID History Tracker

The main objectives of this project are:
1. Retrieve latest coordinates for a list of ENSEMBL IDs
2. Find the complete (GENCODE) annotation history of a gene ID
3. Find the complete (GENCODE) annotation history of a gene symbol
4. Generate putative mappings for deprecated/retired ENSEMBL IDs
5. Build locally an index of GENCODE annotations

The project is being developed in Python 3. A conda environment could be setup as:

```bash
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda create -n genehtrack -c bioconda -c conda-forge -c anaconda nomkl pandas matplotlib networkx bedtools graphviz pygraphviz python=3.7
```

---
## TODO:
- Generate command line script
- Generate full tutorial
- Include functions to download complete Gencode and build the index

---
# Showcase
![map](img/map_example.png)


### Import
```python
from genehtrack import GHTracker

GH = GHTracker(index_path='index')
GH.load_index()
GH.load_graph()
```

### Query by gene list
```python
gene_list = [
    'ENSG00000230544',
    'ENSG00000240453',
    'ENSG00000254100',
    'ENSG00000259484',
    'ENSG00000259758',
    'ENSG00000260940',
    'ENSG00000261496',
    'ENSG00000261685',
    'ENSG00000266584',
    'ENSG00000122952'
]
GH.query_by_geneid_list(gene_list)
```

| idx  | gene_shortid     | gene_id             | gene_name       | hgnc_id     | gene_type             | coord                         | gencode  | ensembl  | assembly  | date        | deprecated | 
|------|------------------|---------------------|-----------------|-------------|-----------------------|-------------------------------|----------|----------|-----------|-------------|------------| 
| 0    | ENSG00000122952  | ENSG00000122952.17  | ZWINT           | HGNC:13195  | protein_coding        | chr10:56357227-56361273(-)    | 37       | 103      | GRCh38    | 2020-12-07  | 0          | 
| 1    | ENSG00000230544  | ENSG00000230544.1   | LINC00453       | NaN         | lincRNA               | chr13:114586640-114588308(+)  | 18       | 73       | GRCh37    | 2013-09-02  | 1          | 
| 2    | ENSG00000240453  | ENSG00000240453.1   | RP11-206L10.10  | NaN         | processed_transcript  | chr1:810109-817712(-)         | 20       | 76       | GRCh38    | 2014-08-26  | 1          | 
| 3    | ENSG00000254100  | ENSG00000254100.1   | AC069120.3      | NaN         | lincRNA               | chr8:38552248-38559020(-)     | 27       | 90       | GRCh38    | 2017-08-01  | 1          | 
| 4    | ENSG00000259484  | ENSG00000259484.1   | RP11-323F24.1   | NaN         | processed_transcript  | chr15:57151866-57210697(-)    | 17       | 72       | GRCh37    | 2013-06-17  | 1          | 
| 5    | ENSG00000259758  | ENSG00000259758.1   | CASC7           | NaN         | lincRNA               | chr8:140520156-140529501(-)   | 21       | 77       | GRCh38    | 2014-09-29  | 1          | 
| 6    | ENSG00000260940  | ENSG00000260940.1   | RP4-575N6.5     | NaN         | sense_overlapping     | chr1:101243158-101243749(+)   | 28       | 92       | GRCh38    | 2018-03-23  | 1          | 
| 7    | ENSG00000261496  | ENSG00000261496.1   | RP13-514E23.1   | NaN         | sense_overlapping     | chr4:86012296-86013874(-)     | 25       | 85       | GRCh38    | 2016-07-15  | 1          | 
| 8    | ENSG00000261685  | ENSG00000261685.2   | RP11-401P9.4    | NaN         | lincRNA               | chr16:50645809-50649249(+)    | 24       | 83       | GRCh38    | 2015-12-03  | 1          | 
| 9    | ENSG00000266584  | ENSG00000266584.1   | AL355149.1      | NaN         | miRNA                 | chr1:16548914-16548987(+)     | 24       | 83       | GRCh38    | 2015-12-03  | 1          | 

