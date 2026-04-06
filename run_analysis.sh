#!/bin/bash

SNP_FILE=$1
MODELS_DIR=$2
OUTPUT_DIR=$3

if [[ -z "$SNP_FILE" ]] || [[ -z "$MODELS_DIR" ]] || [[ -z "$OUTPUT_DIR" ]]; then
    echo "Usage: bash run_analysis.sh <snp_file> <models_dir> <output_dir>"
    echo "Example: bash run_analysis.sh snps.txt /path/to/SMC_models /path/to/results"
    exit 1
fi

LABEL=$(basename "${OUTPUT_DIR}")

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
CHROM_SIZE="/nfs/baldar/quanyiz/genome/hg38/hg38.chrom.sizes"
GENOME_FA="/nfs/baldar/quanyiz/genome/hg38/hg38.fa"
SCORER_PATH="/nfs/baldar/quanyiz/app/variant-scorer/src"

echo "=================================================="
echo "Starting Analysis for: ${LABEL}"
echo "=================================================="

mkdir -p "${OUTPUT_DIR}"

# ================= GPU Detection & Management =================
echo ">>> Detecting Available GPUs"
if ! command -v nvidia-smi &> /dev/null; then
    echo "Error: nvidia-smi not found. Cannot detect GPUs."
    exit 1
fi

FREE_MEM_THRESHOLD=2000 

get_hardware_free_gpus() {
    nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits | \
    awk -v thresh="$FREE_MEM_THRESHOLD" -F', ' '$2 < thresh {print $1}'
}

mapfile -t ALL_GPUS < <(nvidia-smi --query-gpu=index --format=csv,noheader)
NUM_GPUS=${#ALL_GPUS[@]}

if [ "$NUM_GPUS" -eq 0 ]; then
    echo "Error: No GPUs available on this system."
    exit 1
fi
echo "Detected $NUM_GPUS total GPU(s) on the system."
echo "=================================================="

# Initialize GPU PID tracker
declare -A gpu_pids
for gpu in "${ALL_GPUS[@]}"; do
    gpu_pids[$gpu]=""
done

# ================= Process Folds in Parallel =================
for fold in {0..4}; do
    model_path="${MODELS_DIR}/fold_${fold}/models/chrombpnet_nobias.h5"
    
    if [ ! -f "${model_path}" ]; then
        echo "    [Warning] Model not found: ${model_path}. Skipping Fold ${fold}..."
        continue
    fi

    if [ -f "${OUTPUT_DIR}/fold_${fold}.variant_predictions.h5" ]; then
         echo "    > Fold ${fold} results exist. Skipping..."
         continue
    fi

    # Find an available GPU
    assigned_gpu=""
    while [ -z "$assigned_gpu" ]; do
        hardware_free_gpus=($(get_hardware_free_gpus))
        
        for gpu in "${hardware_free_gpus[@]}"; do
            pid=${gpu_pids[$gpu]}
            if [ -z "$pid" ] || ! kill -0 "$pid" 2>/dev/null; then
                assigned_gpu=$gpu
                break
            fi
        done
        
        if [ -z "$assigned_gpu" ]; then
            # Wait for any job to finish or sleep briefly before checking again
            wait -n 2>/dev/null || sleep 10
        fi
    done
    
    echo "  > Assigned [Fold ${fold}] to GPU ${assigned_gpu}"
    
    # Run SHAP and Scoring sequentially on the assigned GPU in the background
    (
        export CUDA_VISIBLE_DEVICES=$assigned_gpu
        
        # Run SHAP
        python ${SCORER_PATH}/variant_shap.py \
            -l ${SNP_FILE} -g ${GENOME_FA} -m ${model_path} -o "${OUTPUT_DIR}/fold_${fold}" \
            -s ${CHROM_SIZE} -sc chrombpnet > /dev/null 2>&1

        # Run Scoring
        python ${SCORER_PATH}/variant_scoring.py \
            -l ${SNP_FILE} -g ${GENOME_FA} -m ${model_path} -o "${OUTPUT_DIR}/fold_${fold}" \
            -s ${CHROM_SIZE} -sc chrombpnet > /dev/null 2>&1
            
        echo "    [Done] Fold ${fold} finished on GPU ${assigned_gpu}"
    ) &
    
    # Track the PID for this GPU
    gpu_pids[$assigned_gpu]=$!
done

# Wait for all background GPU jobs to complete before proceeding
echo "=================================================="
echo "Waiting for all parallel fold computations to finish..."
wait
echo "All fold computations complete."
echo "=================================================="

# ================= Post-Processing =================
echo "  > Running Ensemble & Filtering..."
python "${SCRIPT_DIR}/filter_snps.py" \
    -d "${OUTPUT_DIR}" \
    -c "${LABEL}" \
    -o "${OUTPUT_DIR}" 

echo "  > Generating Plots..."
mkdir -p "${OUTPUT_DIR}/plots"
python "${SCRIPT_DIR}/plot_snps.py" \
    -s "${SNP_FILE}" \
    -d "${OUTPUT_DIR}" \
    -c "${LABEL}" \
    -o "${OUTPUT_DIR}/plots"

echo "=================================================="
echo "Done! Results: ${OUTPUT_DIR}"
echo "=================================================="