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
    norm_x = x - np.max(x, axis=1, keepdims=True)
    exp_x = np.exp(temp * norm_x)
    return exp_x / np.sum(exp_x, axis=1, keepdims=True)

def plot_single_variant(args_tuple):
    row, a1_pred_single, a2_pred_single, a1_shap_single, a2_shap_single, output_dir, cluster_name = args_tuple
    
    variant_id = row['variant_id']
    outfile = os.path.join(output_dir, f"{variant_id}.pdf")

    center_shap = a1_shap_single.shape[0] // 2
    center_pred = a1_pred_single.shape[0] // 2
    flank = 150
    
    try:
        a1_shap = a1_shap_single[center_shap-flank:center_shap+flank]
        a2_shap = a2_shap_single[center_shap-flank:center_shap+flank]
        a1_pred = a1_pred_single[center_pred-flank:center_pred+flank]
        a2_pred = a2_pred_single[center_pred-flank:center_pred+flank]

        if a1_shap.size == 0 or a1_pred.size == 0:
            print(f"Error: {variant_id} 切片为空! 检查数据形状: SHAP={a1_shap_single.shape}")
            return

        v_min = np.nanmin([np.nanmin(a1_shap), np.nanmin(a2_shap)])
        v_max = np.nanmax([np.nanmax(a1_shap), np.nanmax(a2_shap)])
        
        padding = (v_max - v_min) * 0.1 if v_max != v_min else 0.05
        ylim = [v_min - padding, v_max + padding]
        if ylim[0] == ylim[1]: ylim[1] += 0.01

        x = np.arange(-flank, flank)
        c = np.zeros(flank * 2)

        fig = figure(figsize=(20, 9))
        
        # Plot 1: Profile Prediction
        ax1 = fig.add_subplot(311)
        ax1.plot(x, c, color='black')
        ax1.axvline(x=0, color='black', ls='--', linewidth=1)
        ax1.set_xlim(-flank, flank)
        ax1.set_title(f"{variant_id} ({row['ref_allele']}/{row['alt_allele']}) --- {cluster_name}", fontsize=18, weight='bold')
        
        ax1.plot(x, a1_pred, color='royalblue', label=f"Ref ({row['ref_allele']})")
        ax1.plot(x, a2_pred, color='firebrick', label=f"Alt ({row['alt_allele']})")
        ax1.legend(prop={'size': 14}, loc='upper right')

        # Plot 2: Ref SHAP Logo
        ax2 = fig.add_subplot(312)
        df1 = pd.DataFrame(a1_shap, columns=['A','C','G','T'])
        logo1 = logomaker.Logo(df1, ax=ax2)
        logo1.ax.set_ylim(ylim)
        logo1.ax.axvline(x=flank, color='black', ls='--', linewidth=1)
        
        ticks = range(0, (flank * 2) + 1, 50)
        logo1.ax.set_xticks(ticks)
        logo1.ax.set_xticklabels([str(i-flank) for i in ticks])
        ax2.text(0.98, 0.90, "Ref SHAP", transform=ax2.transAxes, size=14, weight='bold', ha='right', bbox=dict(facecolor='white', alpha=0.8))

        # Plot 3: Alt SHAP Logo
        ax3 = fig.add_subplot(313)
        df2 = pd.DataFrame(a2_shap, columns=['A','C','G','T'])
        logo2 = logomaker.Logo(df2, ax=ax3)
        logo2.ax.set_ylim(ylim)
        logo2.ax.axvline(x=flank, color='black', ls='--', linewidth=1)
        
        logo2.ax.set_xticks(ticks)
        logo2.ax.set_xticklabels([str(i-flank) for i in ticks])
        ax3.text(0.98, 0.90, "Alt SHAP", transform=ax3.transAxes, size=14, weight='bold', ha='right', bbox=dict(facecolor='white', alpha=0.8))

        plt.subplots_adjust(hspace=0.4)
        plt.savefig(outfile, format='pdf')
        plt.close(fig)
        
    except Exception as e:
        print(f"Error plotting {variant_id}: {str(e)}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--snp_file", required=True)
    parser.add_argument("-d", "--base_dir", required=True)
    parser.add_argument("-c", "--cluster", required=True)
    parser.add_argument("-o", "--output_dir", required=True)
    args = parser.parse_args()

    shap_inputs = pd.read_table(args.snp_file, names=['chr', 'pos', 'ref_allele', 'alt_allele', 'variant_id'])
    
    print(f"Loading data from {args.base_dir}")
    folds = [f"fold_{i}" for i in range(5)]
    
    preds_a1, preds_a2 = [], []
    shaps_a1, shaps_a2 = [], []
    
    valid_count = 0
    for fold in folds:
        pred_file = os.path.join(args.base_dir, f"{fold}.variant_predictions.h5")
        shap_file = os.path.join(args.base_dir, f"{fold}.variant_shap.counts.h5")
        
        if os.path.exists(pred_file) and os.path.exists(shap_file):
            print(f"  > Reading {fold}...")
            with h5py.File(pred_file, 'r') as f:
                c1 = np.array(f['observed']['allele1_pred_counts']) # shape (N, 1)
                c2 = np.array(f['observed']['allele2_pred_counts']) # shape (N, 1)
                p1 = np.array(f['observed']['allele1_pred_profiles']) # shape (N, 1000)
                p2 = np.array(f['observed']['allele2_pred_profiles']) # shape (N, 1000)
                
                preds_a1.append(c1 * softmax(p1))
                preds_a2.append(c2 * softmax(p2))

            shap_data = dd.io.load(shap_file)
            alleles = np.array(shap_data['alleles'])
            
            seq_shap = np.array(shap_data['projected_shap']['seq'])
            seq_shap = np.transpose(seq_shap, (0, 2, 1))
            
            shaps_a1.append(seq_shap[alleles == 0])
            shaps_a2.append(seq_shap[alleles == 1])
            valid_count += 1

    if valid_count == 0:
        print("Error: No valid H5 files found in the directory.")
        return

    print(f"Ensembling {valid_count} folds...")
    a1_pred_mean = np.mean(preds_a1, axis=0)
    a2_pred_mean = np.mean(preds_a2, axis=0)
    a1_shap_mean = np.mean(shaps_a1, axis=0)
    a2_shap_mean = np.mean(shaps_a2, axis=0)

    os.makedirs(args.output_dir, exist_ok=True)
    
    tasks = []
    for idx in range(min(len(shap_inputs), a1_shap_mean.shape[0])):
        tasks.append((
            shap_inputs.iloc[idx], 
            a1_pred_mean[idx], a2_pred_mean[idx], 
            a1_shap_mean[idx], a2_shap_mean[idx], 
            args.output_dir, args.cluster
        ))
    
    print(f"Generating {len(tasks)} plots using 20 cores...")
    with Pool(20) as p: 
        p.map(plot_single_variant, tasks)
        
    print("Done.")

if __name__ == "__main__":
    main()