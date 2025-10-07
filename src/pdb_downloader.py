"""
PDB Structure Download Module
"""

import requests
from pathlib import Path

def download_pdb_structure(pdb_id, output_path):
    """
    Download structure file from RCSB PDB
    
    Args:
        pdb_id (str): PDB ID
        output_path (Path): Output file path
    """
    url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    
    print(f"Downloading PDB structure: {pdb_id}")
    response = requests.get(url)
    response.raise_for_status()
    
    with open(output_path, 'wb') as f:
        f.write(response.content)
    
    print(f"PDB structure saved: {output_path}")
    return output_path

def validate_pdb_id(pdb_id):
    """Validate PDB ID format"""
    if len(pdb_id) != 4:
        raise ValueError(f"Invalid PDB ID: {pdb_id}, should be 4 characters")
    return pdb_id.upper()