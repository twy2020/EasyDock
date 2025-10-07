#!/bin/bash

echo "🔧 Setting up molecular docking environment..."

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "❌ Conda not found, please install Anaconda or Miniconda first"
    exit 1
fi

# Create conda environment
ENV_NAME="molecular_docking"
if conda env list | grep -q "$ENV_NAME"; then
    echo "✅ Environment '$ENV_NAME' already exists"
else
    echo "Creating conda environment: $ENV_NAME"
    conda create -n $ENV_NAME python=3.9 -y
fi

# Activate environment
eval "$(conda shell.bash hook)"
conda activate $ENV_NAME

# Install Conda packages
echo "Installing Conda packages..."
conda install -c conda-forge openbabel rdkit -y
conda install -c conda-forge mdanalysis -y
conda install -c conda-forge openmm pdbfixer -y

# Install Python dependencies
echo "Installing Python dependencies..."
pip install -r requirements.txt

# Install additional Python packages
pip install py3Dmol prolif

# Download smina
echo "Downloading smina..."
SMINA_PATH="./smina"
if [ ! -f "$SMINA_PATH" ]; then
    wget -q https://sourceforge.net/projects/smina/files/smina.static/download -O $SMINA_PATH
    chmod +x $SMINA_PATH
    echo "✅ smina downloaded and set as executable"
else
    echo "✅ smina already exists"
fi

# Create necessary directories
mkdir -p docking_results

echo "🎉 Environment setup completed!"
echo ""
echo "Use the following command to activate the environment:"
echo "conda activate $ENV_NAME"
echo ""
echo "Run molecular docking like:"
echo "python src/main.py --pdb_id 8WRF --smiles 'Cn1c2cccc(=O)c-2nc2ccccc21'"