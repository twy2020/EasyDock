"""
Visualization Module
Generate 2D and 3D visualization results
"""

import logging
from pathlib import Path
import subprocess

logger = logging.getLogger(__name__)

def generate_2d_interaction(complex_pdb, output_dir):
    """
    Generate 2D interaction diagram
    
    Args:
        complex_pdb: Complex PDB file
        output_dir: Output directory
    """
    try:
        import MDAnalysis as mda
        import prolif as plf
        from prolif.plotting.network import LigNetwork
        from rdkit import Chem
        
        logger.info("Generating 2D interaction diagram...")
        
        # Load complex
        u = mda.Universe(str(complex_pdb))
        
        # Select ligand and protein
        lig = u.select_atoms("resname UNL")
        prot = u.select_atoms("protein")
        
        if len(lig) == 0:
            logger.warning("Ligand (UNL) not found, trying other residue names")
            # Try common ligand residue names
            for resname in ["LIG", "DRG", "INH"]:
                lig = u.select_atoms(f"resname {resname}")
                if len(lig) > 0:
                    logger.info(f"Found ligand: {resname}")
                    break
        
        if len(lig) == 0:
            logger.error("No ligand found, cannot generate 2D interaction diagram")
            return
        
        # Create ligand molecule object
        lmol = plf.Molecule.from_mda(lig)
        
        # Calculate interaction fingerprint
        fp = plf.Fingerprint()
        fp.run(u.trajectory, lig, prot)
        df = fp.to_dataframe(return_atoms=True)
        
        # Generate network diagram
        net = LigNetwork.from_ifp(
            df, lmol,
            kind="aggregate",
            threshold=0.3,
            rotation=270,
        )
        
        # Save image
        output_file = output_dir / "2d_interactions.html"
        net.save(str(output_file))
        logger.info(f"2D interaction diagram saved: {output_file}")
        
    except ImportError as e:
        logger.warning(f"Cannot generate 2D interaction diagram, missing dependency: {e}")
    except Exception as e:
        logger.error(f"Failed to generate 2D interaction diagram: {e}")

def create_pymol_session(complex_pdb, output_dir):
    """
    Create PyMOL session file
    
    Args:
        complex_pdb: Complex PDB file
        output_dir: Output directory
    """
    try:
        # Create PyMOL script
        pml_script = output_dir / "create_session.pml"
        
        script_content = f"""# PyMOL script - automatically generated
load {complex_pdb}, complex

# Select ligand and surrounding residues
select ligand, resn UNL
select surroundings, byres ligand around 5

# Display settings
hide everything
show sticks, ligand or surroundings
hide sticks, elem H

# Coloring
color yellow, ligand and name C*
color cyan, surroundings and name C*

# Label settings
set label_size, 20
set label_color, grey

# Background and style
bg_color white
set cartoon_transparency, 0.5
show cartoon

# View settings
zoom ligand or surroundings

# Save session
save {output_dir / "3d_session.pse"}
"""
        
        with open(pml_script, 'w') as f:
            f.write(script_content)
        
        # Execute PyMOL (if available)
        try:
            cmd = ["pymol", "-c", str(pml_script)]
            subprocess.run(cmd, capture_output=True, check=True)
            logger.info(f"PyMOL session generated: {output_dir / '3d_session.pse'}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.warning("PyMOL not available, only generating script file")
        
        # Cleanup
        pml_script.unlink()
        
    except Exception as e:
        logger.error(f"Failed to create PyMOL session: {e}")

def generate_3d_view(complex_pdb, output_dir):
    """
    Generate interactive 3D view
    
    Args:
        complex_pdb: Complex PDB file
        output_dir: Output directory
    """
    try:
        import py3Dmol
        from rdkit import Chem
        
        logger.info("Generating 3D interactive view...")
        
        # Read complex
        mol = Chem.MolFromPDBFile(str(complex_pdb))
        if mol is None:
            logger.warning("Cannot read PDB file to generate 3D view")
            return
        
        # Convert to PDB string
        pdb_block = Chem.MolToPDBBlock(mol)
        
        # Create view
        view = py3Dmol.view(width=600, height=400)
        view.addModel(pdb_block, 'pdb')
        
        # Set styles
        view.setStyle({'cartoon': {'color': 'spectrum'}})
        view.setStyle({'resn': 'UNL'}, {'stick': {}})
        view.addSurface(py3Dmol.VDW, {'opacity': 0.7, 'color': 'white'}, {'resn': 'UNL'})
        
        view.zoomTo()
        
        # Handle newline issues in PDB block
        escaped_pdb_block = pdb_block.replace('\n', '\\n').replace("'", "\\'")
        
        # Save as HTML - avoid backslash issues in f-string
        html_content = f"""<html>
<head>
    <title>3D Molecular View</title>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
</head>
<body>
    <div style="height: 400px; width: 600px; position: relative;" class='viewer_3Dmoljs' data-pdb='{escaped_pdb_block}' data-backgroundcolor='0xffffff' data-style='cartoon:color=spectrum|stick:resn=UNL' data-surface='opacity:0.7;color:white;resn:UNL' data-labelres='resn:UNL' data-select='resn:UNL' data-spin='on'></div>
</body>
</html>"""
        
        output_file = output_dir / "3d_view.html"
        with open(output_file, 'w') as f:
            f.write(html_content)
        
        logger.info(f"3D interactive view saved: {output_file}")
        
    except ImportError:
        logger.warning("py3Dmol not installed, skipping 3D interactive view")
    except Exception as e:
        logger.error(f"Failed to generate 3D view: {e}")

def create_complex_structure(receptor_pdb, docked_pdb, output_complex):
    """
    Create receptor-ligand complex structure
    
    Args:
        receptor_pdb: Receptor PDB file
        docked_pdb: Docking result PDB file
        output_complex: Output complex file
    """
    try:
        # Simple text merging method
        complex_lines = []
        
        # Read receptor
        with open(receptor_pdb, 'r') as f:
            for line in f:
                if line.startswith(('ATOM', 'HETATM', 'TER')):
                    complex_lines.append(line)
        
        # Read docked ligand (only take first pose)
        with open(docked_pdb, 'r') as f:
            in_pose = False
            pose_count = 0
            
            for line in f:
                if line.startswith('MODEL'):
                    pose_count += 1
                    if pose_count == 1:  # Only take first pose
                        in_pose = True
                        continue
                    else:
                        break
                
                if line.startswith('ENDMDL'):
                    in_pose = False
                    break
                    
                if in_pose and line.startswith(('ATOM', 'HETATM')):
                    complex_lines.append(line)
        
        # Add end marker
        complex_lines.append("END\n")
        
        # Save complex
        with open(output_complex, 'w') as f:
            f.writelines(complex_lines)
        
        logger.info(f"Complex structure saved: {output_complex}")
        
    except Exception as e:
        logger.error(f"Failed to create complex structure: {e}")

def generate_visualizations(receptor_pdb, docked_pdb, output_complex, output_dir):
    """
    Generate all visualization results
    
    Args:
        receptor_pdb: Receptor PDB file
        docked_pdb: Docking result PDB file
        output_complex: Output complex file path
        output_dir: Output directory
    """
    logger.info("Starting to generate visualization results...")
    
    # 1. Create complex structure
    create_complex_structure(receptor_pdb, docked_pdb, output_complex)
    
    # 2. Generate 2D interaction diagram
    generate_2d_interaction(output_complex, output_dir)
    
    # 3. Generate PyMOL session
    create_pymol_session(output_complex, output_dir)
    
    # 4. Generate 3D interactive view
    generate_3d_view(output_complex, output_dir)
    
    # 5. Convert docking results to SDF format
    from docking import convert_docked_to_sdf
    docked_sdf = output_dir / "Docked.sdf"
    convert_docked_to_sdf(docked_pdb, docked_sdf)
    
    logger.info("Visualization results generation completed")