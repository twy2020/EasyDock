"""
Ligand Molecule Preparation Module
Handles SMILES to 3D structure conversion and conformation generation
"""

import subprocess
from pathlib import Path
import logging
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import multiprocessing
import sys

logger = logging.getLogger(__name__)

def smiles_to_3d(smiles, output_sdf, generate_3d=True):
    """
    Convert SMILES to 3D structure
    
    Args:
        smiles: SMILES string
        output_sdf: Output SDF file
        generate_3d: Whether to generate 3D structure
    """
    try:
        # Use RDKit for processing
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string")
        
        # Add hydrogen atoms
        mol = Chem.AddHs(mol)
        
        if generate_3d:
            # Generate 3D coordinates
            AllChem.EmbedMolecule(mol, AllChem.ETKDG())
            # Energy minimization
            AllChem.UFFOptimizeMolecule(mol)
        
        # Save as SDF
        writer = Chem.SDWriter(str(output_sdf))
        writer.write(mol)
        writer.close()
        
        logger.info(f"Ligand 3D structure generated: {output_sdf}")
        
    except Exception as e:
        logger.warning(f"RDKit processing failed: {e}, attempting to use OpenBabel")
        # Fallback option: use OpenBabel
        _smiles_to_3d_obabel(smiles, output_sdf)

def _smiles_to_3d_obabel(smiles, output_sdf):
    """Use OpenBabel to generate 3D structure from SMILES"""
    # Create temporary SMILES file
    temp_smi = output_sdf.with_suffix('.smi')
    with open(temp_smi, 'w') as f:
        f.write(smiles)
    
    # Use OpenBabel for conversion
    cmd = ["obabel", str(temp_smi), "-osdf", "-O", str(output_sdf), "--gen3d"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        raise RuntimeError(f"OpenBabel conversion failed: {result.stderr}")
    
    temp_smi.unlink()  # Delete temporary file
    logger.info(f"Generated ligand 3D structure using OpenBabel: {output_sdf}")

def generate_single_conformations(mol, n, name):
    """Generate multiple conformations for a single molecule - independent function for serialization"""
    if mol is None:
        return None, [], name
        
    mol_with_h = Chem.AddHs(mol)
    try:
        ps = AllChem.ETKDG()
        ps.pruneRmsThresh = 0.5
        ps.numThreads = 0
        
        ids = AllChem.EmbedMultipleConfs(mol_with_h, n, ps)
        for conf_id in ids:
            AllChem.UFFOptimizeMolecule(mol_with_h, confId=conf_id)
        return mol_with_h, list(ids), name
    except Exception as e:
        logger.warning(f"Conformation generation failed {name}: {e}")
        # Generate single conformation as fallback
        try:
            AllChem.EmbedMolecule(mol_with_h)
            AllChem.UFFOptimizeMolecule(mol_with_h)
            return mol_with_h, [0], name
        except:
            return None, [], name

def generate_conformations_single_thread(input_sdf, output_sdf, num_conformers=5):
    """
    Single-threaded generation of multiple conformations (avoid multiprocessing serialization issues)
    
    Args:
        input_sdf: Input SDF file
        output_sdf: Output multi-conformation SDF file
        num_conformers: Number of conformations
    """
    # Read input molecules
    input_mols = [x for x in Chem.SDMolSupplier(str(input_sdf), removeHs=False)]
    if not input_mols:
        raise ValueError("Cannot read input SDF file")
    
    writer = Chem.SDWriter(str(output_sdf))
    
    for i, mol in enumerate(input_mols):
        if mol is None:
            continue
            
        name = f"ligand_{i}"
        logger.info(f"Generating {num_conformers} conformations for molecule {name}...")
        
        try:
            result_mol, conf_ids, _ = generate_single_conformations(mol, num_conformers, name)
            if result_mol and conf_ids:
                result_mol.SetProp("_Name", name)
                for conf_id in conf_ids:
                    writer.write(result_mol, confId=conf_id)
                logger.info(f"Generated {len(conf_ids)} conformations for molecule {name}")
            else:
                logger.warning(f"Conformation generation failed for molecule {name}, using original structure")
                # Write original molecule as fallback
                mol.SetProp("_Name", name)
                writer.write(mol)
                
        except Exception as e:
            logger.error(f"Conformation generation exception for molecule {name}: {e}")
            # Write original molecule as fallback
            mol.SetProp("_Name", name)
            writer.write(mol)
    
    writer.close()
    
    # Verify output
    try:
        output_mols = list(Chem.SDMolSupplier(str(output_sdf)))
        if not output_mols:
            logger.error("Conformation generation failed, no valid output")
            # Create backup: copy input file
            import shutil
            shutil.copy2(input_sdf, output_sdf)
            logger.info(f"Using original structure as backup: {output_sdf}")
        else:
            logger.info(f"Successfully generated {len(output_mols)} conformations: {output_sdf}")
    except Exception as e:
        logger.error(f"Output file verification failed: {e}")
        # Create backup: copy input file
        import shutil
        shutil.copy2(input_sdf, output_sdf)
        logger.info(f"Using original structure as backup: {output_sdf}")

def generate_conformations_multi_process(input_sdf, output_sdf, num_conformers=5, num_threads=None):
    """
    Multi-process generation of multiple conformations (original solution, may fail in some environments)
    """
    if num_threads is None:
        num_threads = max(1, multiprocessing.cpu_count() - 1)
    
    # Read input molecules
    input_mols = [x for x in Chem.SDMolSupplier(str(input_sdf), removeHs=False)]
    if not input_mols:
        raise ValueError("Cannot read input SDF file")
    
    try:
        from concurrent import futures
        
        # Multi-process conformation generation
        with futures.ProcessPoolExecutor(max_workers=num_threads) as executor:
            jobs = []
            for i, mol in enumerate(input_mols):
                if mol:
                    name = f"ligand_{i}"
                    job = executor.submit(generate_single_conformations, mol, num_conformers, name)
                    jobs.append(job)
            
            # Collect results and write to file
            writer = Chem.SDWriter(str(output_sdf))
            for job in futures.as_completed(jobs):
                try:
                    mol, ids, name = job.result()
                    if mol and ids:
                        mol.SetProp("_Name", name)
                        for conf_id in ids:
                            writer.write(mol, confId=conf_id)
                except Exception as e:
                    logger.error(f"Conformation generation task failed: {e}")
            writer.close()
        
        # Verify output
        output_mols = list(Chem.SDMolSupplier(str(output_sdf)))
        if not output_mols:
            raise RuntimeError("Multi-process conformation generation failed, no valid output")
        
        logger.info(f"Multi-process generated {len(output_mols)} conformations: {output_sdf}")
        
    except Exception as e:
        logger.warning(f"Multi-process conformation generation failed: {e}, falling back to single-thread")
        generate_conformations_single_thread(input_sdf, output_sdf, num_conformers)

def generate_conformations(input_sdf, output_sdf, num_conformers=5, use_multiprocessing=True):
    """
    Generate multiple conformations
    
    Args:
        input_sdf: Input SDF file
        output_sdf: Output multi-conformation SDF file
        num_conformers: Number of conformations
        use_multiprocessing: Whether to use multiprocessing
    """
    logger.info(f"Starting conformation generation, target count: {num_conformers}")
    
    if use_multiprocessing:
        try:
            generate_conformations_multi_process(input_sdf, output_sdf, num_conformers)
        except Exception as e:
            logger.warning(f"Multi-processing failed, using single-thread: {e}")
            generate_conformations_single_thread(input_sdf, output_sdf, num_conformers)
    else:
        generate_conformations_single_thread(input_sdf, output_sdf, num_conformers)

def prepare_ligand(smiles, output_sdf, output_conformers, num_conformers=5):
    """
    Complete ligand preparation workflow
    
    Args:
        smiles: SMILES string
        output_sdf: Output single-conformation SDF file
        output_conformers: Output multi-conformation SDF file
        num_conformers: Number of conformations
    """
    logger.info("Starting ligand preparation...")
    
    # 1. SMILES to 3D structure conversion
    smiles_to_3d(smiles, output_sdf)
    
    # 2. Generate multiple conformations (default uses single-thread to avoid serialization issues)
    generate_conformations(output_sdf, output_conformers, num_conformers, use_multiprocessing=False)
    
    # 3. Visualize first conformation (optional)
    try:
        mols = list(Chem.SDMolSupplier(str(output_conformers)))
        if mols:
            img = Draw.MolToImage(mols[0], size=(300, 300))
            img_path = output_conformers.with_suffix('.png')
            img.save(img_path)
            logger.info(f"Ligand structure diagram saved: {img_path}")
    except Exception as e:
        logger.warning(f"Ligand visualization failed: {e}")
    
    logger.info("Ligand preparation completed")