# EasyDock - Automated Molecular Docking Tool
[README-CN](https://github.com/twy2020/EasyDock/blob/main/README-CN.md)

[![Python 3.9+](https://img.shields.io/badge/python-3.9+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)
[![Molecular Docking](https://img.shields.io/badge/domain-Molecular%20Docking-orange.svg)](https://en.wikipedia.org/wiki/Molecular_docking)

## 📖 Project Overview

EasyDock is a localized automated molecular docking tool, modified from the Free_Cloud_Docking (AutoDock colab) project. It has been adapted and optimized for molecular docking box settings, simplifying the complex molecular docking workflow into a single command. Developed based on Smina, it supports full protein coverage docking and multi-conformation search, providing researchers with a convenient and efficient molecular docking solution.

**Original Project Author**: https://github.com/quantaosun  
**Author**: Teng  
**Email**: tenwonyun@gmail.com  
**GitHub**: [https://github.com/twy2020](https://github.com/twy2020?tab=repositories)

## ✨ Core Features

### 🎯 Full Protein Coverage Docking
- **Automatic Docking Box Calculation**: Determines the optimal docking region based on protein structure
- **No Manual Selection Required**: Avoids bias in active site selection common in traditional docking
- **Comprehensive Exploration**: Searches for potential binding sites across the entire protein surface

### 🔬 Multi-Conformation Search
- **Multiple Ligand Conformations**: Supports simultaneous docking of multiple ligand conformations
- **Conformation Generation**: Automatically generates 3D conformations from SMILES
- **Conformation Analysis**: Provides conformation energy distribution and structural diversity analysis

### 🛠️ Automated Workflow
- **One-Click Operation**: From PDB ID and SMILES to complete results
- **Smart Preprocessing**: Automatically repairs protein structures and prepares ligands
- **Format Conversion**: Automatically handles PDB, PDBQT, SDF, and other formats

### 📊 Rich Visualization
- **2D Interaction Diagrams**: Detailed protein-ligand interaction analysis
- **3D Interactive Views**: Web-based 3D molecular viewer
- **PyMOL Sessions**: Professional structural biology analysis session files

![dock](https://gitlab.igem.org/2025/software-tools/yau-china/-/raw/main/EasyDock/pic/dock.png)
![result](https://gitlab.igem.org/2025/software-tools/yau-china/-/raw/main/EasyDock/pic/result.png)

## 🚀 Quick Start

### Environment Requirements
- Linux/macOS/Windows (WSL2 recommended for Windows)
- Python 3.9+
- 4GB+ RAM
- 10GB+ disk space

### One-Click Installation

```bash
# Clone the project
git clone https://github.com/twy2020/EasyDock.git
cd EasyDock

# Run the automatic setup script
chmod +x setup_environment.sh
./setup_environment.sh
```

### Manual Installation

```bash
# Create and activate conda environment
conda create -n easydock python=3.9 -y
conda activate easydock

# Install dependencies
conda install -c conda-forge openbabel rdkit mdanalysis openmm pdbfixer -y
pip install -r requirements.txt

# Download smina
wget https://sourceforge.net/projects/smina/files/smina.static/download -O smina
chmod +x smina
```

## 📖 Usage

### Basic Usage

```bash
# Activate environment
conda activate easydock

# Run docking (with default parameters)
python src/main.py --pdb_id 8WRF --smiles "Cn1c2cccc(=O)c-2nc2ccccc21"
```

### Advanced Parameters

```bash
# Custom docking parameters
python src/main.py \
  --pdb_id 8WRF \
  --smiles "Cn1c2cccc(=O)c-2nc2ccccc21" \
  --ligand_name LIG \
  --work_dir my_docking_results \
  --exhaustiveness 64 \
  --num_modes 200
```

### Parameter Description

| Parameter | Description | Default |
|-----------|-------------|---------|
| `--pdb_id` | PDB database identifier (Required) | - |
| `--smiles` | Ligand SMILES string (Required) | - |
| `--ligand_name` | Residue name of ligand in PDB | LIG |
| `--work_dir` | Working directory | docking_results |
| `--exhaustiveness` | Search intensity (8-128) | 32 |
| `--num_modes` | Number of generated conformations | 100 |
| `--config` | Configuration file path | - |

### Configuration File

Create a `config.yaml` file for batch configuration:

```yaml
# Input parameters
pdb_id: "8WRF"
smiles: "Cn1c2cccc(=O)c-2nc2ccccc21"
ligand_name: "LIG"

# Working directory
work_dir: "docking_results"

# Docking parameters
exhaustiveness: 64
num_modes: 200
energy_range: 5

# Ligand preparation
num_conformations: 10
```

Run with configuration file:
```bash
python src/main.py --config config.yaml
```

## 🔧 Workflow

EasyDock automatically executes the following complete workflow:

### 1. 📥 Data Acquisition
- Download protein structure from RCSB PDB
- Validate PDB ID and file integrity

### 2. 🧬 Protein Preparation
- Separate protein and ligand
- Structure repair and optimization
- Add hydrogen atoms and charges
- Convert to PDBQT format

### 3. ⚗️ Ligand Preparation
- SMILES to 3D structure conversion
- Multi-conformation generation and optimization
- Format conversion and validation

### 4. 🎯 Molecular Docking
- Automatic docking box calculation
- Parallel multi-conformation docking
- Energy scoring and ranking

### 5. 📊 Result Analysis
- Conformation energy analysis
- Interaction analysis
- Multiple format outputs

### 6. 👁️ Visualization
- 2D interaction diagrams
- 3D interactive views
- PyMOL session files

## 📁 Output Files

After completion, the working directory contains files similar to:

```
docking_results/
├── 8WRF-receptor.pdb              # Prepared protein structure
├── receptor.pdbqt                 # Protein for docking
├── small.sdf                      # Ligand 3D structure
├── small_conformation.sdf         # Multi-conformation ligand
├── Dockted.pdb                    # Docking results
├── Dockted.sdf                    # Converted results
├── complex_prepared.pdb           # Complex structure
├── 2d_interactions.html           # 2D interaction diagram
├── 3d_view.html                   # 3D interactive view
├── 3d_session.pse                 # PyMOL session file
└── Docked.log                     # Detailed docking log
```

## 🔍 Result Interpretation

### Docking Scores
- **Binding Energy**: Negative values indicate favorable binding (unit: kcal/mol)
- **RMSD**: Structural differences between conformations
- **Conformation Ranking**: Sorted by binding energy

### Interaction Analysis
- Hydrogen bonds, hydrophobic interactions, π-π stacking, etc.
- Key residue identification
- Binding mode analysis

## 🐛 Troubleshooting

### Common Issues

**Q: Cannot find smina executable**  
A: Ensure smina is downloaded and located in PATH, or use the full path

**Q: RDKit import error**  
A: Install using conda: `conda install -c conda-forge rdkit`

**Q: PDB download failure**  
A: Check network connection and PDB ID validity

**Q: Insufficient memory during docking**  
A: Reduce `num_modes` parameter or use a machine with more memory

### Log Debugging

View detailed logs for error information:
```bash
tail -f docking_results/Docked.log
```

## 📚 Technical Details

### Algorithm Core
- **Docking Engine**: Smina (AutoDock Vina fork)
- **Conformation Generation**: RDKit ETKDG method
- **Protein Processing**: PDBFixer + OpenBabel
- **Visualization**: Py3Dmol + Prolif + MDAnalysis

### Performance Optimization
- Automatic parallel processing
- Memory-efficient data structures
- Incremental result saving

## 🤝 Contributing

We welcome community contributions! Please follow these steps:

1. Fork the project
2. Create a feature branch (`git checkout -b feature/AmazingFeature`)
3. Commit your changes (`git commit -m 'Add some AmazingFeature'`)
4. Push to the branch (`git push origin feature/AmazingFeature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- **Smina Team**: For providing an excellent molecular docking engine
- **RDKit Community**: For the powerful cheminformatics toolkit
- **MDAnalysis**: For molecular dynamics analysis tools
- **PyMOL**: For professional molecular visualization software

## 📞 Support & Contact

If you encounter issues or have suggestions:

- 📧 Email: tenwonyun@gmail.com
- 🐛 [Submit an Issue](https://github.com/twy2020/EasyDock/issues)
- 💬 Discussions: [GitHub Discussions](https://github.com/twy2020/EasyDock/discussions)

## 📊 Citation

If you use EasyDock in your research, please cite:

```bibtex
@software{easydock2024,
  title = {EasyDock: Automated Molecular Docking Pipeline},
  author = {Teng},
  year = {2024},
  url = {https://github.com/twy2020/EasyDock},
  note = {Local automated molecular docking tool with full protein coverage}
}

```
