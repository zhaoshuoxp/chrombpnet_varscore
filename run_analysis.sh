#!/bin/bash

SNP_FILE=$1
CELL_TYPE=$2
MODELS_DIR=$3   
BASE_OUT_DIR=$4 

if [[ -z "$SNP_FILE" ]] || [[ -z "$CELL_TYPE" ]] || [[ -z "$MODELS_DIR" ]] || [[ -z "$BASE_OUT_DIR" ]]; then
    echo "Usage: bash run_analysis.sh <snp_file> <cell_type> <models_dir> <output_dir>"
    echo "Example: bash run_analysis.sh my_snps.txt SMC /path/to/models /path/to/my_results"
    exit 1
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
CHROM_SIZE="/nfs/baldar/quanyiz/genome/hg38/hg38.chrom.sizes"
GENOME_FA="/nfs/baldar/quanyiz/genome/hg38/hg38.fa"
SCORER_PATH="/nfs/baldar/quanyiz/app/variant-scorer/src"

echo "=================================================="
echo "Starting Analysis for ${CELL_TYPE}"
echo "Model Dir:  ${MODELS_DIR}"
echo "Output Dir: ${BASE_OUT_DIR}/${CELL_TYPE}"
echo "=================================================="

FINAL_OUT_DIR="${BASE_OUT_DIR}/${CELL_TYPE}"

for fold in {0..4}; do
    model_path="${MODELS_DIR}/fold_${fold}/models/chrombpnet_nobias.h5"
    fold_out_dir="${FINAL_OUT_DIR}/fold_${fold}"
    
    if [ ! -f "${model_path}" ]; then
        echo "    [Warning] Model not found: ${model_path}. Skipping Fold ${fold}..."
        continue
    fi

    if [ -f "${fold_out_dir}/variant_scores.tsv" ]; then
         echo "    > Fold ${fold} output exists. Skipping..."
         continue
    fi

    echo "    > Processing Fold ${fold}..."
    mkdir -p "${fold_out_dir}"

    python ${SCORER_PATH}/variant_shap.py \
        -l ${SNP_FILE} -g ${GENOME_FA} -m ${model_path} -o ${fold_out_dir} \
        -s ${CHROM_SIZE} -sc chrombpnet > /dev/null 2>&1

    python ${SCORER_PATH}/variant_scoring.py \
        -l ${SNP_FILE} -g ${GENOME_FA} -m ${model_path} -o ${fold_out_dir} \
        -s ${CHROM_SIZE} -sc chrombpnet > /dev/null 2>&1
done

echo "  > Running Ensemble & Filtering..."
python "${SCRIPT_DIR}/filter_snps.py" \
    -d "${FINAL_OUT_DIR}" \
    -c "${CELL_TYPE}" \
    -o "${FINAL_OUT_DIR}" 

echo "  > Generating Plots..."
plot_out_dir="${FINAL_OUT_DIR}/plots"
mkdir -p "${plot_out_dir}"
python "${SCRIPT_DIR}/plot_snps.py" \
    -s "${SNP_FILE}" \
    -d "${FINAL_OUT_DIR}" \
    -c "${CELL_TYPE}" \
    -o "${plot_out_dir}"

echo "=================================================="
echo "All Analysis Finished! Results in: ${FINAL_OUT_DIR}"
echo "=================================================="