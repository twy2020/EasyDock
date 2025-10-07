"""
Receptor Protein Preparation Module
Handles PDB file splitting, repair, and format conversion
"""

import subprocess
import tempfile
from pathlib import Path
import logging

logger = logging.getLogger(__name__)

def split_receptor_ligand(pdb_file, output_receptor, ligand_name="LIG"):
    """
    Split receptor and ligand
    
    Args:
        pdb_file: Input PDB file
        output_receptor: Output receptor file
        ligand_name: Ligand residue name
    """
    try:
        # Use MDAnalysis for splitting (alternative to PyMOL)
        from MDAnalysis import Universe
        import MDAnalysis as mda
        
        u = Universe(pdb_file)
        
        # Select protein as receptor
        protein = u.select_atoms("protein")
        # Select ligand - try multiple possible ligand names
        ligand = None
        possible_ligand_names = [ligand_name, "UNL", "UNK", "DRG", "INH", "LGA"]
        
        for lig_name in possible_ligand_names:
            ligand = u.select_atoms(f"resname {lig_name}")
            if len(ligand) > 0:
                logger.info(f"Found ligand: {lig_name}")
                break
        
        # Save receptor
        if len(protein) > 0:
            protein.write(output_receptor)
            logger.info(f"Receptor saved: {output_receptor}")
        else:
            raise ValueError("No protein structure found")
            
        # Save ligand (if needed)
        if ligand and len(ligand) > 0:
            ligand_file = Path(pdb_file).parent / f"native_{ligand_name}.pdb"
            ligand.write(str(ligand_file))
            logger.info(f"Native ligand saved: {ligand_file}")
            return ligand_file
        else:
            logger.warning(f"No ligand found, will use default docking box location")
            # Create an empty ligand file, automatic detection will be used during docking
            ligand_file = Path(pdb_file).parent / f"native_{ligand_name}.pdb"
            with open(ligand_file, 'w') as f:
                f.write("")  # Create empty file
            return ligand_file
            
    except ImportError:
        logger.warning("MDAnalysis not installed, using simple text processing")
        return _split_receptor_ligand_simple(pdb_file, output_receptor, ligand_name)

def _split_receptor_ligand_simple(pdb_file, output_receptor, ligand_name):
    """Simple text processing method to split receptor and ligand"""
    receptor_lines = []
    ligand_lines = []
    
    with open(pdb_file, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                if ligand_name in line:
                    ligand_lines.append(line)
                else:
                    receptor_lines.append(line)
            else:
                receptor_lines.append(line)
    
    # Save receptor
    with open(output_receptor, 'w') as f:
        f.writelines(receptor_lines)
    
    # Save ligand
    if ligand_lines:
        ligand_file = Path(pdb_file).parent / f"native_{ligand_name}.pdb"
        with open(ligand_file, 'w') as f:
            f.writelines(ligand_lines)
        return ligand_file
    
    return None

def fix_receptor_structure(input_pdb, output_pdb, pH=7.0):
    """
    Use PDBFixer to repair receptor structure
    
    Args:
        input_pdb: Input PDB file
        output_pdb: Output repaired PDB file
        pH: pH value for adding hydrogens
    """
    try:
        from pdbfixer import PDBFixer
        from simtk.openmm.app import PDBFile
        
        logger.info(f"Repairing receptor structure: {input_pdb}")
        
        fixer = PDBFixer(filename=str(input_pdb))
        
        # Repair missing residues
        fixer.findMissingResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        
        # Add missing hydrogen atoms
        fixer.addMissingHydrogens(pH)
        
        # Save repaired structure
        with open(output_pdb, 'w') as f:
            PDBFile.writeFile(fixer.topology, fixer.positions, f)
            
        logger.info(f"Repaired receptor saved: {output_pdb}")
        
    except ImportError:
        logger.warning("PDBFixer not installed, skipping structure repair step")
        # If no PDBFixer, directly copy the file
        import shutil
        shutil.copy2(input_pdb, output_pdb)

def convert_to_pdbqt(input_pdb, output_pdbqt, molecule_type="protein"):
    """
    Use OpenBabel to convert format to PDBQT
    
    Args:
        input_pdb: Input PDB file
        output_pdbqt: Output PDBQT file
        molecule_type: Molecule type (protein/ligand)
    """
    try:
        if molecule_type == "protein":
            # For receptor: remove water molecules, add charges
            cmd = ["obabel", str(input_pdb), "-xr", "-O", str(output_pdbqt)]
        else:
            # For ligand
            cmd = ["obabel", str(input_pdb), "-O", str(output_pdbqt)]
        
        result = subprocess.run(cmd, capture_output=True, text=True)
        
        if result.returncode != 0:
            logger.error(f"OpenBabel conversion failed: {result.stderr}")
            raise RuntimeError(f"PDBQT conversion failed: {result.stderr}")
        
        if not output_pdbqt.exists():
            raise FileNotFoundError(f"Output file not generated: {output_pdbqt}")
            
        logger.info(f"Successfully converted to PDBQT: {output_pdbqt}")
        
    except FileNotFoundError:
        logger.error("OpenBabel not installed, cannot convert format")
        raise

def prepare_receptor(input_pdb, output_receptor, output_pdbqt, ligand_name="LIG"):
    """
    Complete receptor preparation workflow
    
    Args:
        input_pdb: Original PDB file
        output_receptor: Output receptor PDB file
        output_pdbqt: Output receptor PDBQT file
        ligand_name: Ligand residue name
    """
    logger.info("Starting receptor preparation...")
    
    # Create temporary directory
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        
        # 1. Split receptor and ligand
        split_receptor = temp_path / "split_receptor.pdb"
        native_ligand = split_receptor_ligand(input_pdb, split_receptor, ligand_name)
        
        # 2. Repair receptor structure
        fixed_receptor = temp_path / "fixed_receptor.pdb"
        fix_receptor_structure(split_receptor, fixed_receptor)
        
        # 3. Convert to PDBQT format
        convert_to_pdbqt(fixed_receptor, output_pdbqt, "protein")
        
        # 4. Copy final receptor file
        import shutil
        shutil.copy2(fixed_receptor, output_receptor)
        
        # 5. Convert native ligand (if exists)
        if native_ligand and native_ligand.exists():
            native_ligand_pdbqt = output_pdbqt.parent / "native_ligand.pdbqt"
            convert_to_pdbqt(native_ligand, native_ligand_pdbqt, "ligand")
    
    logger.info("Receptor preparation completed")