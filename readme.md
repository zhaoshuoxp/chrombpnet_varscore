- # ChromBPNet Variant Analysis Pipeline

  This pipeline automates the process of variant scoring, statistical filtering, and visualization using **ChromBPNet**. It allows for flexible analysis by manually specifying the model fold directory and the desired output location.

  ## Prerequisites & Environment

  Ensure you have activated the Conda environment containing `chrombpnet`, `deepdish`, `logomaker`, and other dependencies.

  Bash

  ```
  # Installation reference: https://github.com/kundajelab/chrombpnet
  conda activate chrombpnet
  ```

  ## Setup

  1. **Script Files**: Ensure the following three scripts are in the same directory:

     - `run_analysis.sh` (Main control script)
     - `filter_snps.py` (Statistical filtering & ensembling)
     - `plot_snps.py` (Visualization & SHAP logo generation)

  2. **Permissions**: Make the main script executable:

     Bash

     ```
     chmod +x run_analysis.sh
     ```

  3. **Configuration**: Open `run_analysis.sh` and verify that the following global paths match your server environment:

     - `CHROM_SIZE`: Path to `hg38.chrom.sizes`.
     - `GENOME_FA`: Path to `hg38.fa`.
     - `SCORER_PATH`: Path to the `variant-scorer/src` directory (containing `variant_shap.py`).

  ## Usage

  The pipeline now requires **four** arguments to provide maximum flexibility:

  Bash

  ```
  ./run_analysis.sh <snp_file> <cell_type> <models_dir> <output_dir>
  ```

  ### Arguments:

  1. **`snp_file`**: Path to the list of variants to score.
  2. **`cell_type`**: A label for the cell type (used for naming folders and plot titles).
  3. **`models_dir`**: The parent directory containing the 5-fold models.
     - *Structure expected*: `<models_dir>/fold_0/models/chrombpnet_nobias.h5` (up to fold_4).
  4. **`output_dir`**: The root directory where results will be saved.

  ### Example:

  Bash

  ```
  ./run_analysis.sh my_snps.txt SMC /nfs/data/models/SMC_v1 /home/user/project/results
  ```

  ------

  ## Input File Format

  The SNP file must be a tab-separated (TSV) file with **no header**, containing these 5 columns:

  Plaintext

  ```
  chr1    2320766    C    T    rs36096196
  chr1    3409348    C    A    rs2493298
  ```

  - **Columns**: Chromosome, Position (1-based), Ref Allele, Alt Allele, Variant ID.

  ------

  ## Output Structure

  Results are organized within your specified `<output_dir>` under a subfolder named after the `<cell_type>`.

  Plaintext

  ```
  <output_dir>/
  └── <cell_type>/
      ├── fold_0/ ... fold_4/         # Raw scoring & SHAP outputs per fold
      ├── <cell_type>_ensemble_scores.tsv  # Mean/Std scores across all 5 folds
      ├── <cell_type>_sig_hits.tsv         # Filtered significant SNPs
      └── plots/                           # Visualization PDFs
          ├── rs36096196.pdf
          └── ...
  ```

  ------

  ## Key Output Files

  - **`\*_ensemble_scores.tsv`**: Contains averaged metrics across the 5 folds, including `logFC`, `JSD` (Jensen-Shannon Divergence), and the combined score `abs_logfc_x_jsd`.
  - **`\*_sig_hits.tsv`**: A filtered subset of the ensemble scores (typically filtered by p-value and effect size).
  - **`plots/\*.pdf`**: Detailed visualizations for each variant, including:
    1. **Predicted Tracks**: Chromatin accessibility profiles for Ref vs. Alt alleles.
    2. **SHAP Logos**: Sequence importance motifs for both Ref and Alt alleles to visualize motif disruption.

  ------

  ## References

  - **ChromBPNet**: https://github.com/kundajelab/chrombpnet
  - **Variant-scorer**: https://github.com/kundajelab/variant-scorer
  - **Paper**: Nair, S., et al. (2022). *"ChromBPNet: A method to correct biases in chromatin accessibility data."*