#!/bin/bash
#SBATCH --job-name=enformer
#SBATCH --partition=gpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --gres=gpu:1
#SBATCH --time=02:00:00
#SBATCH --output=enformer_upstream/enformer_upstream.out
#SBATCH --error=enformer_upstream/enformer_upstream.err

echo "Running on host: $(hostname)"
echo "Allocated GPU(s): $CUDA_VISIBLE_DEVICES"

# Load environment modules
module purge
module load legacy/CentOS7
module load gcc/11.3.0
module load cuda/11.7.1
module load cudnn/8.4.0.27-11.6
module load conda

# Activate conda env
eval "$(conda shell.bash hook)"
conda activate enformer

# Show which GPU TensorFlow sees
echo "Checking GPU availability..."
python -c "import tensorflow as tf; print('GPUs:', tf.config.list_physical_devices('GPU'))"

# Set thread env vars (optional)
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export TF_NUM_INTRAOP_THREADS=2
export TF_NUM_INTEROP_THREADS=2

start_time=$(date +%s)

# Run Enformer code
python ~/enformer_script.py

end_time=$(date +%s)
runtime=$((end_time - start_time))

echo "Total runtime: ${runtime} seconds"