#!/usr/bin/env python3
"""
Molecular Docking Module
Smina-based molecular docking functionality (supports full protein coverage + multi-conformation search + cross-platform compatibility)
"""

import os
import sys
import logging
import subprocess
import shutil  # For cross-platform executable file search
import numpy as np
from pathlib import Path
from typing import Union, Tuple  # Key: Import Union (type union) and Tuple (tuple type)

# Configure logging (uniformly named "docking" for global log management)
logger = logging.getLogger("docking")
logger.setLevel(logging.INFO)  # Default log level
if not logger.handlers:
    # Add console output handler
    console_handler = logging.StreamHandler()
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    console_handler.setFormatter(formatter)
    logger.addHandler(console_handler)


def find_smina() -> str:
    """
    Cross-platform search for smina executable (supports Linux/macOS/Windows)
    Returns:
        str: Absolute path to smina executable
    Raises:
        RuntimeError: Raised when smina is not found
    """
    # 1. First check common installation paths (covers mainstream systems)
    possible_paths = [
        # Linux common paths
        "/usr/local/bin/smina",
        "/usr/bin/smina",
        # macOS common paths
        "/opt/homebrew/bin/smina",
        "/usr/local/sbin/smina",
        # Windows common paths (requires manual installation and path addition)
        "C:/Program Files/smina/smina.exe",
        "C:/Program Files (x86)/smina/smina.exe",
        # Rely on system PATH if no path specified
        "smina",
        "smina.static"  # Static compiled version (Linux/macOS)
    ]

    # Check priority paths
    for path in possible_paths:
        if os.path.exists(path) and os.access(path, os.X_OK):
            logger.info(f"✅ Found smina at priority path: {path}")
            return path

    # 2. Cross-platform search system PATH (relies on shutil.which, built-in Python 3.3+)
    smina_path = shutil.which("smina")
    if smina_path:
        abs_smina_path = os.path.abspath(smina_path)
        logger.info(f"✅ Found smina in system PATH: {abs_smina_path}")
        return abs_smina_path

    # 3. Throw clear error if not found
    raise RuntimeError(
        "❌ smina executable not found, please install following these steps:\n"
        "Linux: sudo apt update && sudo apt install smina\n"
        "macOS: brew install smina (requires Homebrew installation first)\n"
        "Windows: 1. Download smina (https://github.com/mwojcikowski/smina/releases)\n"
        "         2. Extract to any directory (e.g., C:/Program Files/smina)\n"
        "         3. Add that directory to system environment variable PATH"
    )


def calculate_protein_bounding_box(receptor_pdb: Union[str, Path]) -> Tuple[float, float, float, float, float, float]:
    """
    Calculate protein geometric center and docking box size that wraps the entire protein from receptor PDB file
    Only keeps protein atoms (ATOM lines), excludes ligands/water (HETATM lines), adds 4Å margin to ensure complete coverage

    Args:
        receptor_pdb: Receptor PDB file path (e.g., docking_results/receptor.pdb)
    Returns:
        tuple: (center_x, center_y, center_z, size_x, size_y, size_z)
    Raises:
        RuntimeError: No atomic coordinates extracted or file format error
        FileNotFoundError: PDB file does not exist
    """
    # Standardize path to string
    receptor_pdb_str = str(receptor_pdb)
    logger.info(f"📊 Starting protein docking box calculation: {receptor_pdb_str}")

    # Check if file exists
    if not os.path.exists(receptor_pdb_str):
        raise FileNotFoundError(f"❌ Receptor PDB file does not exist: {receptor_pdb_str}")
    if os.path.getsize(receptor_pdb_str) == 0:
        raise RuntimeError(f"❌ Receptor PDB file is empty: {receptor_pdb_str}")

    # Store protein atomic coordinates (only ATOM lines, exclude HETATM)
    coords = []
    with open(receptor_pdb_str, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            # PDB format specification: ATOM lines represent protein atoms, HETATM represents ligands/water/cofactors
            if line.startswith('ATOM'):
                try:
                    # Extract coordinates from PDB fixed columns (x:30-38, y:38-46, z:46-54)
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                    coords.append([x, y, z])
                except ValueError:
                    # Skip lines with format exceptions (don't interrupt process)
                    logger.warning(f"⚠️  Skipping PDB line {line_num} (coordinate format exception): {line[:50]}...")
                    continue

    # Check if valid coordinates were extracted
    if not coords:
        raise RuntimeError(
            f"❌ No protein atomic coordinates extracted from PDB file\n"
            f"Possible reasons: 1. PDB file only contains HETATM (no protein) 2. File format corrupted 3. Encoding error"
        )

    # Convert to numpy array for calculation
    coords_np = np.array(coords)
    # Geometric center (average of x/y/z)
    center_x = np.mean(coords_np[:, 0])
    center_y = np.mean(coords_np[:, 1])
    center_z = np.mean(coords_np[:, 2])
    # Docking box size (max-min in each dimension + 4Å margin to avoid edge atoms exceeding)
    min_x, max_x = np.min(coords_np[:, 0]), np.max(coords_np[:, 0])
    min_y, max_y = np.min(coords_np[:, 1]), np.max(coords_np[:, 1])
    min_z, max_z = np.min(coords_np[:, 2]), np.max(coords_np[:, 2])
    size_x = (max_x - min_x) + 4.0
    size_y = (max_y - min_y) + 4.0
    size_z = (max_z - min_z) + 4.0

    # Log output key parameters (for verification)
    logger.info(f"✅ Protein geometric center: ({center_x:.2f}, {center_y:.2f}, {center_z:.2f}) Å")
    logger.info(f"✅ Docking box size: {size_x:.2f} × {size_y:.2f} × {size_z:.2f} Å (with 4Å margin)")

    return (center_x, center_y, center_z, size_x, size_y, size_z)


def run_docking(
    receptor_pdb: Union[str, Path],
    receptor_pdbqt: Union[str, Path],
    ligand_sdf: Union[str, Path],
    output_dir: Union[str, Path],
    native_ligand: Union[str, Path, None] = None,  # Note: None should also be included in Union
    exhaustiveness: int = 32,
    num_modes: int = 100,
    energy_range: float = 5.0
) -> str:
    """
    Run molecular docking (core functionality): dynamic docking box + multi-conformation search + result logging

    Args:
        receptor_pdb: Receptor PDB file (for calculating docking box, must contain protein atoms)
        receptor_pdbqt: Receptor PDBQT file (smina docking input format, needs preparation)
        ligand_sdf: Ligand SDF file (contains 3D structure and conformations, needs preparation)
        output_dir: Output directory (automatically created, stores docking results and logs)
        native_ligand: Native ligand PDBQT (optional, for reference position, currently unused)
        exhaustiveness: Search exhaustiveness (32-64 recommended, higher values more thorough but slower)
        num_modes: Number of generated conformations (default 100, covers more possible sites)
        energy_range: Energy range for keeping conformations (kcal/mol, default 5, keeps more candidates)
    Returns:
        str: Absolute path to docking result PDB file (contains num_modes conformations)
    Raises:
        RuntimeError: Docking execution failed or invalid parameters
        FileNotFoundError: Input files do not exist
    """
    logger.info("\n" + "="*60)
    logger.info("🎯 Starting molecular docking (full protein coverage + multi-conformation search)")
    logger.info("="*60)

    # 1. Preprocess input parameters (standardize path format + create output directory)
    # Convert to string to avoid Path object compatibility issues
    receptor_pdb_str = str(receptor_pdb)
    receptor_pdbqt_str = str(receptor_pdbqt)
    ligand_sdf_str = str(ligand_sdf)
    output_dir_str = str(output_dir)
    # Create output directory (if it doesn't exist)
    os.makedirs(output_dir_str, exist_ok=True)

    # 2. Validate input file existence
    input_files = [
        ("Receptor PDB", receptor_pdb_str),
        ("Receptor PDBQT", receptor_pdbqt_str),
        ("Ligand SDF", ligand_sdf_str)
    ]
    for file_desc, file_path in input_files:
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"❌ {file_desc} file does not exist: {file_path}")
        if os.path.getsize(file_path) == 0:
            raise RuntimeError(f"❌ {file_desc} file is empty: {file_path}")

    # 3. Calculate protein docking box (call utility function)
    try:
        center_x, center_y, center_z, size_x, size_y, size_z = calculate_protein_bounding_box(receptor_pdb_str)
    except Exception as e:
        raise RuntimeError(f"❌ Failed to calculate docking box: {str(e)}") from e

    # 4. Find smina executable
    try:
        smina_path = find_smina()
    except Exception as e:
        raise RuntimeError(f"❌ Failed to initialize smina: {str(e)}") from e

    # 5. Define output file paths
    output_pdb = os.path.join(output_dir_str, "Docked.pdb")  # Docking results (multi-conformation)
    output_log = os.path.join(output_dir_str, "Docked.log")  # Docking log (contains energy/parameters)
    logger.info(f"📁 Output result path: {os.path.abspath(output_pdb)}")
    logger.info(f"📁 Output log path: {os.path.abspath(output_log)}")

    # 6. Build smina docking command (core parameters)
    cmd = [
        smina_path,
        "--seed", "0",  # Fixed seed to ensure reproducible results
        "-r", receptor_pdbqt_str,  # Receptor (PDBQT format)
        "-l", ligand_sdf_str,       # Ligand (SDF format)
        "-o", output_pdb,          # Output multi-conformation PDB
        "--log", output_log,        # Output detailed log
        # Search parameters (multi-conformation + high exhaustiveness)
        "--exhaustiveness", str(exhaustiveness),  # Search exhaustiveness (default 32)
        "--num_modes", str(num_modes),            # Number of generated conformations (default 100)
        "--energy_range", str(energy_range),      # Keep conformations with energy difference ≤5
        # Dynamic docking box parameters (calculated from protein)
        "--center_x", f"{center_x:.2f}",
        "--center_y", f"{center_y:.2f}",
        "--center_z", f"{center_z:.2f}",
        "--size_x", f"{size_x:.2f}",
        "--size_y", f"{size_y:.2f}",
        "--size_z", f"{size_z:.2f}"
    ]

    # Print final execution command (for debugging)
    logger.info(f"\n⚙️  Executing docking command: \n{' '.join(cmd)}")

    # 7. Execute docking (capture output + error handling)
    try:
        result = subprocess.run(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            encoding='utf-8',
            timeout=3600  # Timeout 1 hour (avoid infinite hanging)
        )
    except subprocess.TimeoutExpired:
        raise RuntimeError(f"❌ Docking timeout (exceeded 1 hour), possibly due to large protein or strict parameters")
    except Exception as e:
        raise RuntimeError(f"❌ Docking process startup failed: {str(e)}") from e

    # 8. Check docking results (non-zero return code indicates failure)
    if result.returncode != 0:
        # Read log (if exists) to assist troubleshooting
        log_content = ""
        if os.path.exists(output_log):
            with open(output_log, 'r', encoding='utf-8') as f:
                log_content = f.read()[-1000:]  # Only read last 1000 characters
        raise RuntimeError(
            f"❌ smina docking execution failed (return code: {result.returncode})\n"
            f" stderr: {result.stderr[:500]}...\n"
            f" Log tail: {log_content if log_content else 'Log not generated'}"
        )

    # 9. Verify output file existence
    if not os.path.exists(output_pdb) or os.path.getsize(output_pdb) == 0:
        raise RuntimeError(f"❌ Docking result file generation failed: {output_pdb} (empty file or not created)")

    # Completion log
    logger.info(f"\n✅ Docking successful!")
    logger.info(f"✅ Generated {num_modes} docking conformations: {os.path.abspath(output_pdb)}")
    logger.info(f"✅ Docking log saved: {os.path.abspath(output_log)}")
    logger.info("="*60 + "\n")

    return os.path.abspath(output_pdb)


def convert_docked_to_sdf(docked_pdb: Union[str, Path], output_sdf: Union[str, Path]) -> str:
    """
    Convert docking result PDB (multi-conformation) to SDF format (compatible with older RDKit versions, preserves conformation ID and hydrogen atoms)

    Args:
        docked_pdb: Docking result PDB file path (smina output, contains MODEL blocks)
        output_sdf: Output SDF file path
    Returns:
        str: Absolute path to output SDF file
    Raises:
        RuntimeError: Format conversion failed or RDKit not installed
        FileNotFoundError: Input PDB file does not exist
    """
    logger.info(f"\n🔄 Starting docking result format conversion (PDB → SDF)")

    # 1. Preprocess paths and check input
    docked_pdb_str = str(docked_pdb)
    output_sdf_str = str(output_sdf)

    # Check input PDB
    if not os.path.exists(docked_pdb_str):
        raise FileNotFoundError(f"❌ Docking result PDB does not exist: {docked_pdb_str}")
    if os.path.getsize(docked_pdb_str) == 0:
        raise RuntimeError(f"❌ Docking result PDB is empty: {docked_pdb_str}")

    # 2. Import RDKit modules (delayed import to avoid affecting other functions if not installed)
    try:
        from rdkit import Chem
        from rdkit.Chem import SDWriter
    except ImportError:
        raise RuntimeError(
            "❌ RDKit not installed, cannot convert format\n"
            "Install command: conda install -y -c conda-forge rdkit or pip install rdkit-pypi"
        ) from None

    # 3. Read multi-conformation PDB (split by MODEL blocks)
    mols = []
    with open(docked_pdb_str, 'r', encoding='utf-8') as f:
        pdb_content = f.read()

    # Split MODEL blocks (smina output multi-conformation format)
    model_blocks = pdb_content.split('MODEL')[1:]  # Skip content before MODEL
    if not model_blocks:  # Single conformation PDB (no MODEL blocks)
        model_blocks = [pdb_content]
    logger.info(f"📊 Detected {len(model_blocks)} docking conformations")

    # 4. Read conformations one by one (compatible with older RDKit versions, no strictParsing parameter)
    for conf_idx, block in enumerate(model_blocks, 1):
        # Restore complete PDB format for conformation (add MODEL and ENDMDL tags)
        conf_block = f"MODEL {conf_idx}\n{block}"
        if 'ENDMDL' not in conf_block:
            conf_block += "\nENDMDL"

        # Read PDB block (core parameters: keep hydrogens + disable sanitize + disable automatic bonding)
        mol = Chem.MolFromPDBBlock(
            conf_block,
            sanitize=False,        # Skip structure validation (docking results may be non-standard)
            removeHs=False,       # Keep hydrogen atoms (required for subsequent analysis)
            proximityBonding=False# Disable automatic distance-based bonding (avoid corrupting docking results)
        )

        if mol is None:
            logger.warning(f"⚠️  Skipping invalid conformation {conf_idx} (possibly PDB format exception)")
            continue

        # Add ID property to conformation (for subsequent filtering)
        mol.SetProp("Conformer_ID", str(conf_idx))
        # Try to extract docking energy (from PDB remarks)
        if 'REMARK VINA RESULT:' in conf_block:
            for line in conf_block.split('\n'):
                if line.startswith('REMARK VINA RESULT:'):
                    try:
                        energy = float(line.split()[3])
                        mol.SetProp("Docking_Energy_kcal/mol", f"{energy:.2f}")
                    except (IndexError, ValueError):
                        pass
        mols.append(mol)

    # 5. Check number of valid conformations
    if len(mols) == 0:
        raise RuntimeError(f"❌ No valid conformations read from PDB: {docked_pdb_str}")
    logger.info(f"✅ Successfully read {len(mols)} valid conformations")

    # 6. Write SDF file (preserve conformation properties)
    try:
        with SDWriter(output_sdf_str) as writer:
            # Preserve key properties (conformation ID + docking energy)
            writer.SetProps(["Conformer_ID", "Docking_Energy_kcal/mol"])
            for mol in mols:
                writer.write(mol)
    except Exception as e:
        raise RuntimeError(f"❌ Failed to write SDF: {str(e)}") from e

    # 7. Verify output SDF
    if not os.path.exists(output_sdf_str) or os.path.getsize(output_sdf_str) == 0:
        raise RuntimeError(f"❌ SDF file generation failed (empty file): {output_sdf_str}")

    # Completion log
    logger.info(f"✅ Format conversion completed! SDF file: {os.path.abspath(output_sdf_str)}")
    return os.path.abspath(output_sdf_str)


def extract_top_poses(docked_sdf: Union[str, Path], output_dir: Union[str, Path], top_n: int = 3) -> list[str]:
    """
    Extract top N lowest energy conformations from docking result SDF (save as separate PDB files)

    Args:
        docked_sdf: Docking result SDF file path (contains multi-conformation and energy properties)
        output_dir: Output directory (stores Top N conformations)
        top_n: Extract top N conformations (default 3)
    Returns:
        list[str]: List of extracted Top N conformation PDB file paths
    Raises:
        RuntimeError: Extraction failed or no valid conformations
        RDKitImportError: RDKit not installed
    """
    logger.info(f"\n🏆 Extracting top {top_n} lowest energy docking conformations")

    # 1. Preprocess paths
    docked_sdf_str = str(docked_sdf)
    output_dir_str = str(output_dir)
    os.makedirs(output_dir_str, exist_ok=True)

    # 2. Import RDKit
    try:
        from rdkit import Chem
    except ImportError:
        raise RuntimeError(
            "❌ RDKit not installed, cannot extract conformations\n"
            "Install command: conda install -y -c conda-forge rdkit or pip install rdkit-pypi"
        ) from None

    # 3. Read SDF (disable sanitize to avoid read failures)
    suppl = Chem.SDMolSupplier(docked_sdf_str, sanitize=False)
    mols = [mol for mol in suppl if mol is not None]
    if not mols:
        raise RuntimeError(f"❌ No valid conformations read from SDF: {docked_sdf_str}")
    logger.info(f"📊 Total detected valid conformations: {len(mols)}")

    # 4. Sort by docking energy (prefer Docking_Energy_kcal/mol property)
    def get_energy(mol: Chem.Mol) -> float:
        """Get molecule docking energy (for sorting)"""
        if mol.HasProp("Docking_Energy_kcal/mol"):
            return float(mol.GetProp("Docking_Energy_kcal/mol"))
        # Compatible with other energy property names (e.g., smina default output)
        elif mol.HasProp("minimizedAffinity"):
            return float(mol.GetProp("minimizedAffinity"))
        elif mol.HasProp("r_i_docking_score"):
            return float(mol.GetProp("r_i_docking_score"))
        else:
            # Sort by conformation ID if no energy property
            return int(mol.GetProp("Conformer_ID")) if mol.HasProp("Conformer_ID") else float('inf')

    # Sort by energy ascending (lower energy is better)
    sorted_mols = sorted(mols, key=get_energy)
    top_mols = sorted_mols[:min(top_n, len(sorted_mols))]

    # 5. Save Top N conformations as PDB
    top_poses = []
    for rank, mol in enumerate(top_mols, 1):
        # Generate output filename (includes rank and energy)
        energy = get_energy(mol)
        output_pdb = os.path.join(output_dir_str, f"Top_{rank}_Energy_{energy:.2f}.pdb")
        output_pdb_abs = os.path.abspath(output_pdb)

        # Save PDB (use with statement to automatically close file, avoid resource leaks)
        try:
            with Chem.PDBWriter(output_pdb_abs) as writer:
                writer.write(mol)
        except Exception as e:
            logger.warning(f"⚠️  Failed to save Top {rank} conformation: {str(e)}")
            continue

        top_poses.append(output_pdb_abs)
        logger.info(f"✅ Top {rank} conformation: {output_pdb_abs} (Energy: {energy:.2f} kcal/mol)")

    # Check extraction results
    if not top_poses:
        raise RuntimeError(f"❌ Failed to extract any conformations (possibly save failure)")

    return top_poses


def analyze_docking_results(docked_sdf: Union[str, Path]) -> dict:
    """
    Analyze docking result SDF: statistics on conformation count, energy distribution, RMSD, etc.

    Args:
        docked_sdf: Docking result SDF file path
    Returns:
        dict: Analysis result dictionary (contains total conformations, energy distribution, RMSD, etc.)
    Raises:
        RuntimeError: Analysis failed or no valid conformations
    """
    logger.info(f"\n📈 Starting docking result analysis")

    # 1. Import RDKit
    try:
        from rdkit import Chem
    except ImportError:
        return {
            "status": "error",
            "message": "RDKit not installed, cannot perform detailed analysis",
            "install_cmd": "conda install -y -c conda-forge rdkit or pip install rdkit-pypi"
        }

    # 2. Read SDF
    docked_sdf_str = str(docked_sdf)
    if not os.path.exists(docked_sdf_str):
        return {"status": "error", "message": f"SDF file does not exist: {docked_sdf_str}"}

    suppl = Chem.SDMolSupplier(docked_sdf_str, sanitize=False)
    mols = [mol for mol in suppl if mol is not None]
    if not mols:
        return {"status": "error", "message": f"No valid conformations read: {docked_sdf_str}"}

    # 3. Extract conformation information
    pose_info_list = []
    energies = []
    for mol in mols:
        pose_info = {
            "conformer_id": mol.GetProp("Conformer_ID") if mol.HasProp("Conformer_ID") else "unknown",
            "energy_kcal/mol": None,
            "rmsd_lb": None,  # RMSD lower bound (smina output)
            "atom_count": mol.GetNumAtoms()
        }

        # Extract energy
        if mol.HasProp("Docking_Energy_kcal/mol"):
            pose_info["energy_kcal/mol"] = round(float(mol.GetProp("Docking_Energy_kcal/mol")), 2)
            energies.append(pose_info["energy_kcal/mol"])
        elif mol.HasProp("minimizedAffinity"):
            pose_info["energy_kcal/mol"] = round(float(mol.GetProp("minimizedAffinity")), 2)
            energies.append(pose_info["energy_kcal/mol"])

        # Extract RMSD (smina output lower bound)
        if mol.HasProp("rmsd_lb"):
            pose_info["rmsd_lb"] = round(float(mol.GetProp("rmsd_lb")), 3)

        pose_info_list.append(pose_info)

    # 4. Calculate energy statistics
    energy_stats = {}
    if energies:
        energy_stats = {
            "min_energy": round(min(energies), 2),
            "max_energy": round(max(energies), 2),
            "avg_energy": round(sum(energies)/len(energies), 2),
            "energy_std": round(np.std(energies), 2) if len(energies) > 1 else 0.0
        }

    # 5. Organize analysis results
    analysis_result = {
        "status": "success",
        "total_conformers": len(mols),
        "energy_statistics": energy_stats,
        "conformers": pose_info_list[:10],  # Only return detailed information for first 10 conformations
        "note": f"Only showing first 10 conformations, complete results contain {len(mols)}"
    }

    # Print key statistical information
    logger.info(f"✅ Analysis completed!")
    logger.info(f"📊 Total conformations: {len(mols)}")
    if energy_stats:
        logger.info(f"📊 Energy distribution: Min {energy_stats['min_energy']} | Avg {energy_stats['avg_energy']} | Max {energy_stats['max_energy']} kcal/mol")

    return analysis_result


# Test code (executed when this file is run directly)
if __name__ == "__main__":
    """
    Test command: python docking.py
    Requires preparing test files in advance (test_data directory):
    - test_data/receptor.pdb (protein PDB)
    - test_data/receptor.pdbqt (protein PDBQT)
    - test_data/ligand.sdf (ligand SDF, contains 3D conformation)
    """
    # Configure test paths
    test_data_dir = Path("test_data")
    test_output_dir = Path("test_docking_output")

    # Check test files
    required_files = [
        test_data_dir / "receptor.pdb",
        test_data_dir / "receptor.pdbqt",
        test_data_dir / "ligand.sdf"
    ]
    missing_files = [str(f) for f in required_files if not f.exists()]
    if missing_files:
        print(f"❌ Missing test files: {missing_files}")
        print(f"⚠️  Please place receptor PDB/PDBQT and ligand SDF in {test_data_dir} directory")
        sys.exit(1)

    try:
        # 1. Run docking
        docked_pdb = run_docking(
            receptor_pdb=test_data_dir / "receptor.pdb",
            receptor_pdbqt=test_data_dir / "receptor.pdbqt",
            ligand_sdf=test_data_dir / "ligand.sdf",
            output_dir=test_output_dir,
            exhaustiveness=16,  # Lower exhaustiveness during testing for speed
            num_modes=10        # Reduce conformation count during testing
        )

        # 2. Convert to SDF
        docked_sdf = convert_docked_to_sdf(
            docked_pdb=docked_pdb,
            output_sdf=test_output_dir / "Docked.sdf"
        )

        # 3. Extract Top 2 conformations
        top_poses = extract_top_poses(
            docked_sdf=docked_sdf,
            output_dir=test_output_dir / "top_poses",
            top_n=2
        )

        # 4. Analyze results
        analysis = analyze_docking_results(docked_sdf=docked_sdf)
        print(f"\n📊 Docking result analysis: {analysis}")

        print(f"\n🎉 All test steps completed! Results saved in: {test_output_dir.absolute()}")

    except Exception as e:
        print(f"\n❌ Test failed: {str(e)}")
        logger.exception("Test failure detailed stack trace")
        sys.exit(1)