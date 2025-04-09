# ClinGen HPO analysis

- Analysis of gene-disease associations based on HPO terms assigned to probands within the ClinGen database
- Uses ClinGen API to query proband HPO terms
- Includes multiple analytic workflows using those HPO terms

### Important packages

Python

- `pyprojroot` - Enables project-based workflow in python, setting root path to root of repository
- `hpo3` - Module for HPO ontology-based calculations. Allows parallelization and multiple similarity metrics
- Graph based modules
  - `networkx`
  - `igraph`
  - `leidenalg`

R

- `here` - Enables project-based workflow in R, setting root path to root of repository
- `ontologySimilarity` - Library with tools for ontology similarity calculations. Non-parallelized and only `lin` and `resnik` based similarity

### Scripts dir

- `scripts/gcep.py` - Contains `gcep` class with functionality to query ClinGen API and generate HPO tables for each proband
- `scripts/gcep_config.py` - Must be created. Should contain the following information:
- `scripts/gcep_scrape.py` - Contains `gcep_scrape` class with functionality to query ClinGen website directly without API. Allows query of all probands without a unique key necessary for each GCEP 
- `scripts/gcep_scrape_pipeline.py` - Contains script to query ClinGen database from command line using an HGNC gene ID. 
- `scripts/create_hpo_distance_object_scrape.py` - Takes compiled output from scraping approach and finds distances between all probands based on HPO sets
- `scripts/process_iuis_table.R` - Processes raw IUIS excel file downloaded from https://iuis.org/committees/iei/ into `data/raw_data`

```python
gcep_url="[ClinGen API URL]"
api_key_pird = "[PIRD GCEP API Key]"
affiliation_pird= "[PIRD GCEP Affiliation ID]"
api_key_scid = "[SCID GCEP API Key]"
affiliation_scid = "[SCID GCEP Affiliation ID]"
```

- `create_hpo_distance_object.py` - Makes a query using a `gcep` object and returns two distance matrices: (1) the distance between probands based on HPO-term sets and (2) the distance between all individual HPO-terms present in the data. Packages these matrices and associated metadata into an `hdf5` file that is saved to `data/`
- `scripts/hpo_openai_embedding.py` - Obtains the OpenAI embedding values for the name of all HPO terms in the HPO database and writes it to the data dir. 

### Data dir

- `hp_[date].obo` - HPO ontology file
- `genes_to_phenotype_[date].txt` - All gene-hpo term associations contained within the HPO database
- `data/clingen_scrape` - Contains results of scraping approach. Run in two iterations, a primary one which most genes and a second which was used to complete any genes that failed after debugging the original script. 
- `data/clingen_scrape/gcep_key.csv` - Compiled csv of all proband data from above 
- `data/clingen_scrape/gcep_key.csv`- Pulls data out of `.pkl` files in same directories which creates a key of which genes are associated with which GCEP committees
- `data/clingen_scrape/clingen_scrape_hpo.h5` - H5 object which contains distance matrix of probands based on HPO sets and also contains associated metadata
- `data/iuis_table.csv` - Table of IUIS genes, there groups, and subgroups
