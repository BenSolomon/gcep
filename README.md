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

```python
gcep_url="[ClinGen API URL]"
api_key_pird = "[PIRD GCEP API Key]"
affiliation_pird= "[PIRD GCEP Affiliation ID]"
api_key_scid = "[SCID GCEP API Key]"
affiliation_scid = "[SCID GCEP Affiliation ID]"
```

- `create_hpo_distance_object.py` - Makes a query using a `gcep` object and returns two distance matrices: (1) the distance between probands based on HPO-term sets and (2) the distance between all individual HPO-terms present in the data. Packages these matrices and associated metadata into an `hdf5` file that is saved to `data/`

### Data dir

- `hp_[date].obo` - HPO ontology file
- `genes_to_phenotype_[date].txt` - All gene-hpo term associations contained within the HPO database
