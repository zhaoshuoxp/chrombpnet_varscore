#!/bin/bash

SNP_FILE=$1
CELL_TYPE=$2

if [[ -z "$SNP_FILE" ]] || [[ -z "$CELL_TYPE" ]]; then
    echo "Usage: bash run_analysis.sh <snp_file> <cell_type>"
    echo "Example: bash run_analysis.sh my_snps.txt SMC"
    exit 1
fi

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
echo "Script directory is: ${SCRIPT_DIR}"

CHROM_SIZE="/nfs/baldar/quanyiz/genome/hg38/hg38.chrom.sizes"
GENOME_FA="/nfs/baldar/quanyiz/genome/hg38/hg38.fa"
SCORER_PATH="/nfs/baldar/quanyiz/app/variant-scorer/src"

MULTIOME_MODEL_ROOT="/nfs/baldar/quanyiz/Multiome/chrombpnet/model"
CZI_MODEL_ROOT="/nfs/baldar/quanyiz/CZI/ATAC/chromBPnet/model"


run_pipeline() {
    local dataset_name=$1  
    local model_root=$2
    local folder_name=$3  
    local output_root="./${dataset_name}_Results"
    
    echo "=================================================="
    echo "Starting ${dataset_name} Analysis for ${CELL_TYPE} (Model Folder: ${folder_name})"
    echo "=================================================="

    local base_out_dir="${output_root}/${CELL_TYPE}" #

    for fold in {0..4}; do
        local model_path="${model_root}/${folder_name}/fold_${fold}/models/chrombpnet_nobias.h5"
        local fold_out_dir="${base_out_dir}/fold_${fold}"
        
        if [ ! -f "${model_path}" ]; then
            echo "    [Warning] Model not found: ${model_path}. Skipping Fold ${fold}..."
            continue
        fi

        if [ -f "${fold_out_dir}/variant_scores.tsv" ]; then
             echo "    > Fold ${fold} output exists. Skipping..."
             continue
        fi

        echo "    > Processing Fold ${fold}..."

        python ${SCORER_PATH}/variant_shap.py \
            -l ${SNP_FILE} -g ${GENOME_FA} -m ${model_path} -o ${fold_out_dir} \
            -s ${CHROM_SIZE} -sc chrombpnet > /dev/null 2>&1

        python ${SCORER_PATH}/variant_scoring.py \
            -l ${SNP_FILE} -g ${GENOME_FA} -m ${model_path} -o ${fold_out_dir} \
            -s ${CHROM_SIZE} -sc chrombpnet > /dev/null 2>&1
    done

    echo "  > Running Ensemble & Filtering..."
    python "${SCRIPT_DIR}/filter_snps.py" \
        -d "${base_out_dir}" \
        -c "${CELL_TYPE}" \
        -o "${base_out_dir}" 

    echo "  > Generating Plots..."
    local plot_out_dir="${base_out_dir}/plots"
    python "${SCRIPT_DIR}/plot_snps.py" \
        -s "${SNP_FILE}" \
        -d "${base_out_dir}" \
        -c "${CELL_TYPE}" \
        -o "${plot_out_dir}"

    echo "  > ${dataset_name} Done. Results in ${base_out_dir}"
}


run_pipeline "Multiome" "${MULTIOME_MODEL_ROOT}" "${CELL_TYPE}"

czi_folder=""

case "${CELL_TYPE}" in
    "SMC")
        czi_folder="Smooth_Muscle_Medial"
        ;;
    "modSMC")
        czi_folder="Smooth_Muscle_Intimal"
        ;;
    "Microvasculature")
        czi_folder="Microvasculature"
        ;;
    "Endothelial"|"Fibroblast"|"Lymphocyte"|"Macrophage")
        czi_folder="${CELL_TYPE}"
        ;;
    *)
        echo "Warning: No CZI mapping found for ${CELL_TYPE}. Skipping CZI run."
        czi_folder=""
        ;;
esac

if [[ -n "$czi_folder" ]]; then
    run_pipeline "CZI" "${CZI_MODEL_ROOT}" "${czi_folder}"
fi

echo "=================================================="
echo "All Analysis Finished!"
echo "=================================================="