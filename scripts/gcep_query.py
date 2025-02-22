import requests
import gcep_config
import pandas as pd
import re


class gcep:
    """
    Class that takes parameters query to the GCEP API and parses response data
    Includes methods to return data in a variety of formats
    """    
    def __init__(self, api_key, gcep_url, status, affiliation, start, end):  
        self.api_key = api_key
        self.gcep_url = gcep_url
        self.params = {
            "target": "gci",
            "status": status,
            "affiliation": affiliation,
            "start": start,
            "end": end
        }
        self.response = self._api_get()
        self.json = self.response.json()
        self.table = pd.DataFrame(self.json)
        self.genes = [x['Gene'] for x in self.json]
        self.n_genes = len(self.genes)
        # self.hpo_table = self._hpo_table()
    
    def _api_get(self):
        """
        Internal method used by __init__ to make a get request to the GCEP API

        Returns:
            _type_: _description_
        """        
        response = requests.get(
            f"{self.gcep_url}/snapshots",
            headers={"x-api-key": self.api_key},
            params=self.params)
        return response    
    
    def _single_gene_proband_table(self, row):
        """
        Internal method used as pandas apply function to format proband data for a single gene
        Parses the proband sub-json format and joins it with the original Gene and Disease data

        Args:
            row (_type_): _description_

        Returns:
            _type_: _description_
        """        
        df = pd.DataFrame(row['probands'])
        gene = row['Gene']
        disease = row['Disease']
        df['Gene'] = gene
        df['Disease'] = disease
        df = df[['Gene', 'Disease'] + [col for col in df.columns if col not in ['Gene', 'Disease']]]
        return df
    
    
    def _format_hpo_string(self, hpo_str):
        """
        Internal method used to reformat how HPO terms are stored in GCEP json data
        Intake format from GCEP JSON follows pattern "HPO term (HPO_ID)"
        Output format is a dictionary with keys "HPO_ID" and "HPO_term"
        This allows the HPO_ID to be used with HPO3 module for further analysis
        
        Args:
            hpo_str (_type_): _description_

        Returns:
            _type_: _description_
        """        
        try: 
            match = re.match(r"(.+?)\s*\((HP:\d+)\)", hpo_str)
            if match:
                hpo_dict = {
                    "HPO_ID": match.group(2),
                    "HPO_term": match.group(1)
                }
                return hpo_dict
            else:
                return hpo_str, None
        except:
            return hpo_str, None
        
    
    def proband_table(self):
        """
        Takes data from self.table and returns a new table with each proband in the dataset.
        Each entry in self.table['probands'] is a json tree of proband data that is converted to a dataframe
        The original Gene and Disease from the original row is added to this new dataframe
        The resulting proband dataframes are concatenated into a single dataframe
        """        
        df = self.table.apply(self._single_gene_proband_table, axis=1)
        df = pd.concat(df.to_list(), axis = 0, ignore_index=True)
        return df
    
    def hpo_table(self):
        """
        Takes data from self.table and returns a new table for each HPO term for every 
        Gene-Disease-proband combination
        Includes application of _format_hpo_string to reformat HPO terms

        Returns:
            _type_: _description_
        """        
        df = self.proband_table()
        df = df[['Gene', 'Disease', 'label','HPO terms']].explode("HPO terms")
        df['HPO terms'] = df['HPO terms'].apply(self._format_hpo_string) 
        df = df.join(pd.json_normalize(df.pop('HPO terms')))
        df = df.dropna(subset=['HPO_ID'])
        return df
                
gcep_query = gcep(
    api_key = gcep_config.api_key_pird, 
    gcep_url = gcep_config.gcep_url,
    status = "approved",
    affiliation=gcep_config.affiliation_pird,
    start = "2024-12-01", 
    end = "2025-01-31"
    )

print(gcep_query.genes)
print(gcep_query.table)
print(gcep_query.proband_table())
print(gcep_query.hpo_table())