import pandas as pd
import numpy as np
import os
import argparse
import sys

def mean_filter_sig_snp(base_dir, cluster_name, output_dir, logfc=0.25, p=0.05):
    os.makedirs(output_dir, exist_ok=True)
    ensemble_file = os.path.join(output_dir, f"{cluster_name}_ensemble_scores.tsv")
    sig_hits_file = os.path.join(output_dir, f"{cluster_name}_sig_hits.tsv")

    NUMERIC_COLS = [
        'allele1_pred_counts', 'allele2_pred_counts',
        'logfc', 'abs_logfc', 
        'jsd', 'original_jsd', 
        'logfc_x_jsd', 'abs_logfc_x_jsd',
        'logfc.pval', 'abs_logfc.pval', 'jsd.pval' 
    ]
    
    dfs = []
    meta_df = None  

    for fold in range(5):
        file_name = f"fold_{fold}.variant_scores.tsv"
        file_path = os.path.join(base_dir, file_name)
        
        if os.path.exists(file_path):
            print(f"  Loading Fold {fold}...")
            try:
                df = pd.read_csv(file_path, sep="\t")
                if 'variant_id' in df.columns:
                    df = df.set_index('variant_id')
                else:
                    print(f"    Error: 'variant_id' column not found in fold {fold}")
                    continue
                    
                cols_to_use = [c for c in NUMERIC_COLS if c in df.columns]
                dfs.append(df[cols_to_use])
                
                if meta_df is None:
                    meta_cols = [c for c in df.columns if c not in NUMERIC_COLS]
                    meta_df = df[meta_cols]
                    
            except Exception as e:
                print(f"    Error reading fold {fold}: {e}")
        else:
            print(f"  Warning: File not found: {file_path}")

    if not dfs:
        print(f"No data found for {cluster_name}. Check directory: {base_dir}")
        return

    print("Calculating ensemble means...")
    ensemble_df = pd.concat(dfs).groupby(level=0).mean()
    
    if meta_df is not None:
        final_df = meta_df.join(ensemble_df, how='inner')
    else:
        final_df = ensemble_df 

    std_df = pd.concat(dfs).groupby(level=0)['logfc'].std()
    final_df['logfc_std'] = std_df

    final_df = final_df.sort_values('abs_logfc_x_jsd', ascending=False)
    final_df.to_csv(ensemble_file, sep="\t")
    print(f"Ensemble results saved to: {ensemble_file}")

    if 'abs_logfc_x_jsd.pval' in final_df.columns and 'abs_logfc' in final_df.columns:
        cond_pval = final_df['abs_logfc_x_jsd.pval'] < p
        cond_logfc = final_df['abs_logfc'] > logfc
        
        sig_df = final_df[cond_pval & cond_logfc].copy()
        sig_df = sig_df.sort_values('abs_logfc_x_jsd', ascending=False)

        percentage = (len(sig_df)/len(final_df))*100 if len(final_df) > 0 else 0
        print(f"Significant hits found: {len(sig_df)} ({percentage:.2f}%)")
        
        sig_df.to_csv(sig_hits_file, sep="\t")
        print(f"Top hits saved to: {sig_hits_file}")
    else:
        print("Required columns for filtering not found.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--base_dir", required=True, help="Directory containing fold folders")
    parser.add_argument("-c", "--cluster", required=True, help="Cluster name")
    parser.add_argument("-o", "--output", required=True, help="Output directory")
    args = parser.parse_args()
    
    mean_filter_sig_snp(args.base_dir, args.cluster, args.output)