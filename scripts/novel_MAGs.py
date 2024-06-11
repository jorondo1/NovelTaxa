#!/usr/bin/env python/3.11.5
# coding=utf-8

import pandas as pd
import argparse
import os 
os.environ['PYTHONUNBUFFERED'] = '1' # print all messages

def main():
    parser = argparse.ArgumentParser(description="Identify novel species-level MAGs from skANI output.")

    # Add arguments
    parser.add_argument('-a', '--ANI_results', type=str, required=True, help='File containing skANI results.')
    parser.add_argument('-m', '--MAG_list', type=str, required=True, help='File with with full path to MAGs (one per lines)', default='tmp/MAG_list.txt')
    parser.add_argument('-c', '--checkm', type=str, required=True, help='Checkm quality report.')
    parser.add_argument('-o', '--outdir', type=str, required=True, help='output directory')
    parser.add_argument('-g', '--gtdb', type=str, required=True, help='GTDB metadata file')
    args = parser.parse_args() # Parse

    checkm50, ani95, ani = identifyNovel(args.ANI_results, args.checkm, args.MAG_list, args.outdir)
    myMAGsBringAllTheBoysToTheYardDamnRightItsBetterThanYours(checkm50, ani95, ani, args.gtdb, args.outdir)

def identifyNovel(ani_file, checkm_file, MAG_file, out):
 # MAGs, quality score and ANI
    MAGs = pd.read_csv(MAG_file, sep="\t", header = None)
    checkm = pd.read_csv(checkm_file, sep='\t')                 
    ani = pd.read_csv(ani_file, sep="\t")[['Ref_file','Query_file', 'ANI']]
    print(f"{len(MAGs)} MAGs evaluated.")
    print(f"{len(ani['Query_file'].drop_duplicates())} MAGs have an ANI >= 80% with at least one GTDB genome.")

    # List MAGs with <95% ANI_results
    ani95 = ani[ani['ANI']>=95.00]['Query_file'].drop_duplicates()
    print(f"{len(ani95)} MAGs share a species cluster (ANI >= 95%) with at least one GTDB genome.")

    checkm['QS'] = checkm['Completeness'] - 5 * checkm['Contamination']

    # List MAGs with QS > 50
    checkm50 = checkm.query('QS>=50.00')[['Name','QS']]
    print(f"{len(checkm50)} MAGs have quality score (QS) over 0.50.")

    exclude = (set(checkm50['Name'].tolist()) - set(ani95.tolist()))
    print(f'{len(exclude)} MAGs are potentially novel species-level MAGs with QS > 0.50.')

    # Filter paths 
    nMAGs = MAGs[MAGs[0].apply(lambda x: any(s in x for s in exclude))]
    nMAGs.to_csv(f"{out}/nMAG_list.txt", index=False, header = False)
    
    return checkm50, ani95, ani
    
def myMAGsBringAllTheBoysToTheYardDamnRightItsBetterThanYours(checkm50, ani95, ani, gtdb_file, out):
    gtdb=pd.read_csv(gtdb_file, sep='\t')
    
    #  Find the Best MAGs with ANI >=95%
    intersect=set(checkm50['Name']).intersection(set(ani95))
    idx = ani[ani['Query_file'].isin(intersect)].groupby('Query_file')['ANI'].idxmax() # extract indexes 

    redunSpec = (
    	ani.loc[idx] # keep intersect
    	.merge(checkm50, left_on='Query_file', right_on='Name', how='left') # add QS
    	.drop(columns=['Name']) # remove dup
    )
    # remove first 3 letters from accession name (they all start with either GB_ or RS_)
    gtdb['accession'] = gtdb['accession'].str[3:]
    # compute QS
    gtdb['QS_ref'] = gtdb['checkm2_completeness'] - 5 * gtdb['checkm2_contamination']

    # Merge both and check best value 
    compRef = redunSpec.merge(gtdb[['accession','QS_ref']], left_on='Ref_file', right_on='accession', how='left').drop(columns='Ref_file')
    compRef.to_csv(f"{out}/betterMAGs.txt", index=False, header = True, sep='\t')
    betterMAGs = compRef[compRef['QS']>compRef['QS_ref']].copy()
    betterMAGs['increase'] = (betterMAGs['QS'] - betterMAGs['QS_ref'])/betterMAGs['QS_ref']
    meanInc = betterMAGs['increase'].mean() * 100
    sdInc = betterMAGs['increase'].std() * 100
    print(f"found {sum(compRef['QS'] > compRef['QS_ref'])} MAG with higher quality score (mean increase {meanInc:.1f} ± {sdInc:.1f}%)")
    print(betterMAGs)

if __name__ == "__main__": # Ensure function call only when script is run directly, not loaded as a module.
    main()
