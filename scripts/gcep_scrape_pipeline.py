import sys
from pyprojroot.here import here
import argparse
import pandas as pd
import gzip

here()
sys.path.insert(0, here('scripts'))
from gcep_scrape import gcep_scrape


def main():
    parser = argparse.ArgumentParser(description="Scrape ClinGen probands")
    parser.add_argument('HGNC_ID', type=str, nargs='?', default="HGNC:12731", help='The gene HGNC ID to scrape')
    parser.add_argument('SAVE_DIR', type=str, nargs='?', default = here(), help='Directory to save output files')
    args = parser.parse_args()
    
    print(f'### STARTING {args.HGNC_ID}', file=sys.stderr)
    
    
    gcep_query = gcep_scrape(args.HGNC_ID)
    print(f'{args.HGNC_ID}\t{gcep_query.valid_entry}')

        
    if gcep_query.valid_entry:
        with gzip.open(f'{args.SAVE_DIR}/{args.HGNC_ID.replace(":", "_")}.pkl.gz', 'wb') as f:
            gcep_query.table.to_pickle(f)
        df_hpo=gcep_query.hpo_table()
        if not df_hpo.empty:
            df_hpo.to_csv(f'{args.SAVE_DIR}/{args.HGNC_ID.replace(':', '_')}_hpo.csv' , index=False)


if __name__ ==  "__main__":
    main()