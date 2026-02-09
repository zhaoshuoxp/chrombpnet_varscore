# ChromBPNet Variant Analysis Pipeline

This pipeline automates the process of variant scoring, statistical filtering, and visualization using **ChromBPNet**. It supports analysis across multiple datasets (e.g., Multiome & CZI) and handles cell-type name mapping automatically.

## Prerequisites & Environment

Before running the pipeline, ensure you have activated the correct Conda environment containing `chrombpnet`, `deepdish`, `logomaker`, etc.

Bash

```
# installation:https://github.com/kundajelab/chrombpnet?tab=readme-ov-file#installation
conda activate chrombpnet
```

## Setup

1. Ensure the following three scripts are in the same directory:

   - `run_analysis.sh` (Main control script)
   - `filter_snps.py` (Statistical filtering script)
   - `plot_snps.py` (Visualization script)

2. Make the main script executable:

   Bash

   ```
   chmod +x run_analysis.sh
   ```

3. **IMPORTANT Configuration**: Open `run_analysis.sh` and ensure the following paths point to the correct locations on your server:

   - `CHROM_SIZE`
   - `GENOME_FA`
   - `SCORER_PATH` (Path to `variant_shap.py` / `variant_scoring.py`)
   - `MULTIOME_MODEL_ROOT` & `CZI_MODEL_ROOT` (Paths to trained models)

## Usage

Run the pipeline using the `run_analysis.sh` script. It requires two arguments: the SNP list file and the Cell Type name.

Bash

```
./run_analysis.sh <snp_file> <cell_type>
```

### Examples

To analyze SMC (Smooth Muscle Cells):

Bash

```
./run_analysis.sh CAD_loci_SNP.txt SMC
```

To analyze modSMC (Modulated SMC):

Bash

```
./run_analysis.sh CAD_loci_SNP.txt modSMC
```

*Note: The script automatically handles mapping between Multiome names (e.g., `SMC`) and CZI names (e.g., `Smooth_Muscle_Medial`).*

## Input File Format

The SNP file must be a tab-separated or space-separated file with **no header**, containing the following 5 columns:

Plaintext

```
chr1    2320766    C    T    rs36096196
chr1    3409348    C    A    rs2493298
chr1    3409946    G    A    rs116710059
```

*Columns: Chromosome, Position (1-based), Ref Allele, Alt Allele, Variant ID.*

## Output Structure

Results are organized by dataset (`Multiome_Results` or `CZI_Results`) and Cell Type.

Plaintext

```
./Multiome_Results/
    └── SMC/
        ├── fold_0/ ... fold_4/         # Raw scoring output
        ├── SMC_ensemble_scores.tsv     # Aggregated scores (Mean/Std across folds)
        ├── SMC_sig_hits.tsv            # Filtered significant SNPs (Top hits)
        └── plots/                      # Visualization PDFs
            ├── rs36096196.pdf
            └── ...
```

### Key Output Files

- **`\*_ensemble_scores.tsv`**: Contains metrics like `logfc`, `jsd`, and the combined score `abs_logfc_x_jsd` averaged across 5 folds.
- **`\*_sig_hits.tsv`**: A subset of SNPs filtered by significance (default: p<0.05, logfc>0.25).
- **`plots/\*.pdf`**: Visualization containing:
  1. Predicted chromatin accessibility profile (Ref vs Alt).
  2. Ref allele SHAP motif logo.
  3. Alt allele SHAP motif logo.

## References

This pipeline utilizes **ChromBPNet** for deep learning-based chromatin accessibility prediction.

- **Chrombpnet Repository**: https://github.com/kundajelab/chrombpnet
- **Variant-scorer Repository**: https://github.com/kundajelab/variant-scorer
- **Paper**: *Nair, S., et al. (2022). "ChromBPNet: A method to correct biases in chromatin accessibility data."*