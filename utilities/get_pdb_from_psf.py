# import libraries
import sys, os
import MDAnalysis as mda

print("Export pdb...")

INPUT_PATH = os.getcwd()
OVERWRITE = 'YES'  #'NO' to skip existing files

# get directories
folders = [f for f in os.listdir(INPUT_PATH) if os.path.isdir(os.path.join(INPUT_PATH, f))]
print(f'[INFO] Found folders: {folders}')
print(f'[INFO] # number of system: {len(folders)}')


for receptor in folders:
    psf = os.path.join(INPUT_PATH, receptor, "build", "structure.psf")
    pdb = os.path.join(INPUT_PATH, receptor, "build", "structure.pdb")

    if not os.path.exists(psf) or not os.path.exists(pdb):
        print(f'[WARNING] File not found: {file_path}')
        if not os.path.exists(psf):
            print(f'Missing: {psf}')
        if not os.path.exists(pdb):
            print(f'Missing: {pdb}')
        continue
    print(f'[INFO] Found both PSF and PDB for {receptor}')
    
    # Extract and save only the protein

    file_out = os.path.join(INPUT_PATH, receptor, f"{receptor}_only_prot.pdb")
    if not os.path.exists(file_out):
        print(f"Saving file: {file_out}")
        u = mda.Universe(psf, pdb)
        protein = u.select_atoms("protein")
        protein.write(file_out, reindex=False)
    elif OVERWRITE == 'YES':
        print(f"Overwriting existing file: {file_out}")
        u = mda.Universe(psf, pdb)
        protein = u.select_atoms("protein")
        protein.write(file_out, reindex=False)
    else:
        print(f"[INFO] Skipping {receptor}: output already exists!")







