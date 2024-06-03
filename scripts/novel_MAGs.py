#!/usr/bin/env python3
# coding=utf-8

import pandas as pd
import argparse
import os 
os.environ['PYTHONUNBUFFERED'] = '1' # print all messages

def main():
    parser = argparse.ArgumentParser(description="Identify novel species-level MAGs from skANI output.")

    # Add arguments
    parser.add_argument('-a', '--ANI_results', type=str, required=True, help='File containing skANI results.', default='output/ANI_results.txt')
    parser.add_argument('-m', '--MAG_list', type=str, help='File with with full path to MAGs (one per lines)', default='tmp/MAG_list.txt')
    parser.add_argument('-c', '--checkm', type=str, help='Checkm quality report.', default='checkm2/quality_report.tsv')

    # Parse the arguments
    args = parser.parse_args()

    # Read MAG list
    MAGs = pd.read_csv(args.MAG_list, sep="\t", header = None)
    print(f"{len(MAGs)} MAGs evaluated.")

    # Import ANI results
    ani = pd.read_csv(args.ANI_results, sep="\t")[['Ref_file','Query_file', 'ANI', 'Align_fraction_query']]
    print(f"{len(ani['Query_file'].drop_duplicates())} MAGs have an ANI >= 80% with at least one GTDB genome.")

    # List MAGs with <95% ANI
    ani95 = ani[ani['ANI']>=95.00]['Query_file'].drop_duplicates().tolist()
    print(f"{len(ani95)} MAGs share a species cluster (ANI >= 95%) with at least one GTDB genome.")

    # Import checkM results
    checkm = pd.read_csv(args.checkm, sep='\t')
    checkm['QS'] = checkm['Completeness'] - 5 * checkm['Contamination']

    # List MAGs with QS > 50
    checkm50 = checkm.query('QS>=50.00')['Name'].tolist()
    print(f"{len(checkm50)} MAGs have quality score (QS) over 0.50.")

    intersect = (set(checkm50) - set(ani95))
    print(f'{len(intersect)} MAGs are potentially novel species-level MAGs with QS > 0.50.')

    # Filter paths 
    nMAGs = MAGs[MAGs[0].apply(lambda x: any(s in x for s in intersect))]
    nMAGs.to_csv('tmp/nMAG_list.txt', index=False, header = False)

# Ensure function call only when script is run directly, not loaded as a module.
if __name__ == "__main__":
    main()
