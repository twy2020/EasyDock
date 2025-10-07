#!/usr/bin/env python3
"""
Molecular Docking Main Program - Local Version
Supports full protein coverage docking and multi-conformation search
"""

import argparse
import os
import sys
from pathlib import Path
import yaml
import logging

# Add src directory to Python path
sys.path.insert(0, str(Path(__file__).parent))

try:
    from pdb_downloader import download_pdb_structure
    from receptor_prep import prepare_receptor
    from ligand_prep import prepare_ligand
    from docking import run_docking  # Use modified docking function
    from visualization import generate_visualizations
except ImportError as e:
    logging.error(f"Module import failed: {e}")
    logging.error("Please ensure all dependencies are correctly installed")
    sys.exit(1)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)


class MolecularDockingPipeline:
    def __init__(self, config):
        self.config = config
        self.work_dir = Path(config.get('work_dir', 'docking_results'))
        self.work_dir.mkdir(exist_ok=True)
        
        # Core parameters (obtained from config, supports command line/config file settings)
        self.pdb_id = config['pdb_id']
        self.smiles = config['smiles']
        self.ligand_name = config.get('ligand_name', 'LIG')
        self.exhaustiveness = config.get('exhaustiveness', 32)  # Search exhaustiveness (default 32)
        self.num_modes = config.get('num_modes', 100)  # Number of docking conformations (default 100)
        self.energy_range = config.get('energy_range', 5)  # Energy range (default 5)
        
        # File path management
        self.file_paths = {
            'raw_pdb': self.work_dir / f"{self.pdb_id}.pdb",
            'receptor_pdb': self.work_dir / f"{self.pdb_id}-receptor.pdb",  # For calculating docking box
            'receptor_pdbqt': self.work_dir / "receptor.pdbqt",  # Receptor for docking
            'native_ligand_pdbqt': self.work_dir / "native_ligand.pdbqt",
            'ligand_sdf': self.work_dir / "small.sdf",
            'ligand_conformations': self.work_dir / "small_conformation.sdf",  # Ligand conformations
            'docked_results': self.work_dir / "Docked.pdb",  # Docking results
            'docked_sdf': self.work_dir / "Docked.sdf",
            'complex_structure': self.work_dir / "complex_prepared.pdb"
        }
    
    def run_pipeline(self):
        """Run complete molecular docking pipeline"""
        print("🚀 Starting molecular docking pipeline...")
        
        try:
            # 1. Download PDB structure
            print("\n📥 Step 1: Downloading PDB structure...")
            self.download_pdb()
            
            # 2. Prepare receptor (generate protein structure with hydrogens)
            print("\n🔧 Step 2: Preparing receptor...")
            self.prepare_receptor()
            
            # 3. Prepare ligand (generate 3D conformations)
            print("\n⚗️ Step 3: Preparing ligand...")
            self.prepare_ligand()
            
            # 4. Execute docking (core modification: pass receptor PDB and multi-conformation parameters)
            print("\n🎯 Step 4: Executing molecular docking...")
            self.run_docking()
            
            # 5. Generate visualization results
            print("\n📊 Step 5: Generating visualization results...")
            self.generate_visualization()
            
            print(f"\n✅ Molecular docking completed! Results saved in: {self.work_dir}")
            
        except Exception as e:
            print(f"❌ Pipeline execution failed: {e}")
            logger.exception("Detailed information of pipeline execution failure")
            sys.exit(1)
    
    def download_pdb(self):
        """Download PDB structure (if it doesn't exist)"""
        if not self.file_paths['raw_pdb'].exists():
            download_pdb_structure(self.pdb_id, self.file_paths['raw_pdb'])
        else:
            print(f"PDB file already exists: {self.file_paths['raw_pdb']}")
    
    def prepare_receptor(self):
        """Prepare receptor protein (ensure generation of PDB file for calculating docking box)"""
        prepare_receptor(
            input_pdb=self.file_paths['raw_pdb'],
            output_receptor=self.file_paths['receptor_pdb'],  # Corrected parameter name
            output_pdbqt=self.file_paths['receptor_pdbqt'],
            ligand_name=self.ligand_name
        )
        # Verify receptor file generation
        if not self.file_paths['receptor_pdb'].exists():
            raise FileNotFoundError(f"Receptor PDB file generation failed: {self.file_paths['receptor_pdb']}")
    
    def prepare_ligand(self):
        """Prepare ligand molecule (generate 3D conformations)"""
        num_confs = self.config.get('num_conformations', 5)  # Ligand's own conformation count
        prepare_ligand(
            smiles=self.smiles,
            output_sdf=self.file_paths['ligand_sdf'],
            output_conformers=self.file_paths['ligand_conformations'],  # Corrected parameter name
            num_conformers=num_confs  # Corrected parameter name
        )
    
    def run_docking(self):
        """Execute molecular docking (core: pass receptor PDB for calculating docking box)"""
        # Call modified run_docking function, pass necessary parameters
        docked_result_str = run_docking(
            receptor_pdb=self.file_paths['receptor_pdb'],  # For calculating docking box
            receptor_pdbqt=self.file_paths['receptor_pdbqt'],  # Receptor for docking
            ligand_sdf=self.file_paths['ligand_conformations'],  # Ligand conformations
            output_dir=self.work_dir,
            native_ligand=self.file_paths['native_ligand_pdbqt'],
            exhaustiveness=self.exhaustiveness,  # Search exhaustiveness
            num_modes=self.num_modes,  # Generate 100 conformations
            energy_range=self.energy_range  # Energy range
        )
        # Convert returned string to Path object
        self.file_paths['docked_results'] = Path(docked_result_str)
    
    def generate_visualization(self):
        """Generate visualization results (2D interaction diagrams and 3D structures)"""
        generate_visualizations(
            receptor_pdb=self.file_paths['receptor_pdb'],
            docked_pdb=self.file_paths['docked_results'],
            output_complex=self.file_paths['complex_structure'],  # Corrected parameter name
            output_dir=self.work_dir
            # Remove redundant generate_2d and generate_3d parameters
        )


def load_config(config_file=None):
    """Load configuration file (default values + user configuration)"""
    # Default configuration (includes docking optimization parameters)
    default_config = {
        'work_dir': 'docking_results',
        'exhaustiveness': 32,         # Search exhaustiveness (recommended 32-64)
        'num_modes': 100,             # Number of docking conformations
        'energy_range': 5,            # Energy range for keeping conformations
        'num_conformations': 5,       # Ligand's own conformation count
        'generate_2d': True,          # Whether to generate 2D interaction diagrams
        'generate_3d': True           # Whether to generate 3D interactive views
    }
    
    # Merge user configuration (if any)
    if config_file and Path(config_file).exists():
        with open(config_file, 'r') as f:
            user_config = yaml.safe_load(f) or {}
        default_config.update(user_config)
    
    return default_config


def main():
    parser = argparse.ArgumentParser(description='Molecular Docking Local Version (Supports full protein coverage + multi-conformation search)')
    parser.add_argument('--pdb_id', required=True, help='PDB ID (e.g.: 8WRF)')
    parser.add_argument('--smiles', required=True, help='Ligand SMILES string')
    parser.add_argument('--ligand_name', default='LIG', help='Ligand three-letter code (used to remove original ligand from PDB)')
    parser.add_argument('--work_dir', default='docking_results', help='Working directory')
    parser.add_argument('--config', help='Configuration file path (YAML format)')
    # New: Allow direct setting of core docking parameters via command line
    parser.add_argument('--exhaustiveness', type=int, help=f'Search exhaustiveness (default 32)')
    parser.add_argument('--num_modes', type=int, help=f'Number of docking conformations (default 100)')
    
    args = parser.parse_args()
    
    # Build configuration (command line arguments have higher priority than config file)
    config = load_config(args.config)
    # Basic configuration
    config.update({
        'pdb_id': args.pdb_id,
        'smiles': args.smiles,
        'ligand_name': args.ligand_name,
        'work_dir': args.work_dir
    })
    # Docking parameters (override if specified via command line)
    if args.exhaustiveness is not None:
        config['exhaustiveness'] = args.exhaustiveness
    if args.num_modes is not None:
        config['num_modes'] = args.num_modes
    
    # Print key configuration (for user confirmation)
    logger.info(f"Key docking parameters: exhaustiveness={config['exhaustiveness']}, num_modes={config['num_modes']}")
    
    # Run pipeline
    pipeline = MolecularDockingPipeline(config)
    pipeline.run_pipeline()


if __name__ == "__main__":
    main()