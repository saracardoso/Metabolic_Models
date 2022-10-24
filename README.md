# Metabolic_Models

With the aid of the scRNAseq data from this [colorectal cancer (CRC) atlas](https://github.com/saracardoso/CRC_ATLAS), we reconstructed, for each sample, cell-type specific models for several T-cell subtypes.

This repository stores all code developed to reconstruct and analyse the models, as well as the models reconstructed.

| **Cell Type**     | **Number of Models** | **State Distribution**        | **CMS Distribution**                      |
|-------------------|----------------------|-------------------------------|-------------------------------------------|
| Cytotoxic CD8     | 16                   | Tumour: 10 Normal Matched: 6  | CMS1: 2 CMS2: 6 CMS3: 2 CMS4: 0 Mixed: 0  |
| Follicular CD4    | 13                   | Tumour: 13 Normal Matched: 0  | CMS1: 3 CMS2: 6 CMS3: 3 CMS4: 0 Mixed: 1  |
| IL17+ CD4         | 13                   | Tumour: 8 Normal Matched: 5   | CMS1: 2 CMS2: 4 CMS3: 1 CMS4: 0 Mixed: 1  |
| Memory CD4        | 33                   | Tumour: 20 Normal Matched: 13 | CMS1: 6 CMS2: 6 CMS3: 3 CMS4: 4 Mixed: 1  |
| Memory CD8        | 30                   | Tumour: 18 Normal Matched: 12 | CMS1: 6 CMS2: 6 CMS3: 3 CMS4: 2 Mixed: 1  |
| Naive CD4         | 30                   | Tumour: 19 Normal Matched: 11 | CMS1: 6 CMS2: 6 CMS3: 2 CMS4: 4 Mixed: 1  |
| Naive CD8         | 11                   | Tumour: 6 Normal Matched: 5   | CMS1: 0 CMS2: 4 CMS3: 2 CMS4: 0 Mixed: 0  |
| Proliferative CD4 | 12                   | Tumour: 12 Normal Matched: 0  | CMS1: 3 CMS2: 4 CMS3: 2 CMS4: 2 Mixed: 1  |
| Proliferative CD8 | 9                    | Tumour: 9 Normal Matched: 0   | CMS1: 3 CMS2: 3 CMS3: 2 CMS4: 0 Mixed: 1  |
| Regulatory CD4    | 29                   | Tumour: 20 Normal Matched: 9  | CMS1: 6 CMS2: 6 CMS3: 3 CMS4: 4 Mixed: 1  |


## Where can the models be found?

The T-cells models can be found, in the *SBML* format, in the following folder of this repository: `0MODELS/CRC_atlas`. The folder is organized as follows:

` /<cell_type>/`

`     /NormalMatched/ OR /Tumour/`

`          /<indiv><sample><cell_type>.xml.gz`


## Structure of the project



<!-- ## How to reference -->



