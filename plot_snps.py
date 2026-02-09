import pandas as pd
import numpy as np
import os
import h5py
import deepdish as dd
import argparse
from matplotlib import pyplot as plt
from matplotlib.pyplot import figure
import logomaker
from multiprocessing import Pool

import matplotlib
matplotlib.use('Agg')

def softmax(x, temp=1):
    norm_x = x - np.mean(x, axis=1, keepdims=True)
    return np.exp(temp*norm_x)/np.sum(np.exp(temp*norm_x), axis=1, keepdims=True)

def plot_single_variant(args_tuple):
    row, allele1_pred_mean, allele2_pred_mean, allele1_shap_mean, allele2_shap_mean, output_dir, cluster_name = args_tuple
    
    variant_id = row['variant_id']
    outfile = os.path.join(output_dir, f"{variant_id}.pdf")

    center_shap = 1057 
    center_pred = 500 
    flank = 150
    
    try:
        a1_shap = allele1_shap_mean[center_shap-flank:center_shap+flank]
        a2_shap = allele2_shap_mean[center_shap-flank:center_shap+flank]
        a1_pred = allele1_pred_mean[center_pred-flank:center_pred+flank]
        a2_pred = allele2_pred_mean[center_pred-flank:center_pred+flank]

        ylim_min = min(a1_shap.min(), a2_shap.min()) * 1.1
        ylim_max = max(a1_shap.max(), a2_shap.max()) * 1.1
        
        if ylim_max == ylim_min: ylim_max += 0.01

        ylim = [ylim_min, ylim_max]
        
        x = np.arange(-flank, flank)
        c = np.zeros(flank * 2)

        fig = figure(figsize=(20, 9))
        
        # Plot 1: Profile Prediction
        ax1 = fig.add_subplot(311)
        plt.plot(x, c, color='black')
        plt.axvline(x=0, color='black', ls='--', linewidth=1)
        ax1.margins(x=0, y=0.1)
        ax1.set_xlim(-flank, flank + 1)
        plt.title(f"{variant_id} ({row['ref_allele']}/{row['alt_allele']}) --- {cluster_name}", fontsize=18, weight='bold')
        
        plt.plot(x, a1_pred, color='royalblue', label=f"ref ({row['ref_allele']})")
        plt.plot(x, a2_pred, color='firebrick', label=f"alt ({row['alt_allele']})")
        plt.legend(prop={'size': 18}, loc='upper right')

        # Plot 2: Ref SHAP Logo
        ax2 = fig.add_subplot(312)
        logo1 = logomaker.Logo(pd.DataFrame(a1_shap, columns=['A','C','G','T']), ax=ax2)
        logo1.ax.set_xlim(0, (flank * 2))
        logo1.ax.set_ylim(ylim)
        ax2.axvline(x=flank, color='black', ls='--', linewidth=1)
     
        ticks = range(0, (flank * 2) + 1, 50)
        logo1.ax.set_xticks(ticks)
        logo1.ax.set_xticklabels([str(i-flank) for i in ticks])
        plt.text(0.98, 0.90, "Ref", transform=ax2.transAxes, size=18, bbox=dict(facecolor='white', alpha=0.5))

        # Plot 3: Alt SHAP Logo
        ax3 = fig.add_subplot(313)
        logo2 = logomaker.Logo(pd.DataFrame(a2_shap, columns=['A','C','G','T']), ax=ax3)
        logo2.ax.set_xlim(0, (flank * 2))
        logo2.ax.set_ylim(ylim)
        ax3.axvline(x=flank, color='black', ls='--', linewidth=1)
        logo2.ax.set_xticks(ticks)
        logo2.ax.set_xticklabels([str(i-flank) for i in ticks])
        plt.text(0.98, 0.90, "Alt", transform=ax3.transAxes, size=18, bbox=dict(facecolor='white', alpha=0.5))

        plt.subplots_adjust(hspace=0.4)
        plt.savefig(outfile, format='pdf')
        plt.close()
        
    except Exception as e:
        print(f"Error plotting {variant_id}: {e}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--snp_file", required=True)
    parser.add_argument("-d", "--base_dir", required=True, help="Directory containing fold_0, fold_1 etc.")
    parser.add_argument("-c", "--cluster", required=True)
    parser.add_argument("-o", "--output_dir", required=True)
    args = parser.parse_args()

    shap_inputs = pd.read_table(args.snp_file, names=['chr', 'pos', 'ref_allele', 'alt_allele', 'variant_id'])
    
    print(f"Loading data for cluster: {args.cluster} from {args.base_dir}")
    folds = [f"fold_{i}" for i in range(5)]
    
    preds_accum = {'a1': [], 'a2': []}
    shap_accum = {'a1': [], 'a2': []}
    
    valid_folds = 0
    for fold in folds:
        pred_file = os.path.join(args.base_dir, f"{fold}.variant_predictions.h5")
        shap_file = os.path.join(args.base_dir, f"{fold}.variant_shap.counts.h5")
        
        if os.path.exists(pred_file) and os.path.exists(shap_file):
            valid_folds += 1
            with h5py.File(pred_file, 'r') as f:
                c1 = np.array(f['observed']['allele1_pred_counts'])
                c2 = np.array(f['observed']['allele2_pred_counts'])
                p1 = np.array(f['observed']['allele1_pred_profiles'])
                p2 = np.array(f['observed']['allele2_pred_profiles'])
                
                # Calculate profile * counts
                preds_accum['a1'].append(c1 * softmax(p1))
                preds_accum['a2'].append(c2 * softmax(p2))

            # Load SHAP
            shap_data = dd.io.load(shap_file)
            alleles = np.array(shap_data['alleles']) # 0=ref, 1=alt
            seq_shap = np.array(shap_data['projected_shap']['seq'])
            
            shap_accum['a1'].append(seq_shap[alleles==0])
            shap_accum['a2'].append(seq_shap[alleles==1])
    
    if valid_folds == 0:
        print("No valid fold files found.")
        return

    print("Calculating means across folds...")
    a1_pred_mean = np.mean(np.array(preds_accum['a1']), axis=0)
    a2_pred_mean = np.mean(np.array(preds_accum['a2']), axis=0)
    a1_shap_mean = np.mean(np.array(shap_accum['a1']), axis=0)
    a2_shap_mean = np.mean(np.array(shap_accum['a2']), axis=0)
    
    print(f"Plotting {len(shap_inputs)} variants...")
    os.makedirs(args.output_dir, exist_ok=True)
    
    tasks = []
    for idx, row in shap_inputs.iterrows():
        tasks.append((
            row, 
            a1_pred_mean[idx], a2_pred_mean[idx], 
            a1_shap_mean[idx], a2_shap_mean[idx], 
            args.output_dir, args.cluster
        ))
    
    with Pool(20) as p: 
        p.map(plot_single_variant, tasks)
        
    print("Plotting done.")

if __name__ == "__main__":
    main()