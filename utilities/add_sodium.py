#!/usr/bin/env python3
"""
Add a sodium ion (Na+) between the carboxylate groups of two acidic residues.

Usage:
    python add_na_between_acids.py input.pdb RES1 NUM1 RES2 NUM2 [OFFSET]

Example:
    python add_na_between_acids.py protein.pdb ASP 70 ASP 111 0.5
"""

import sys
import math
import numpy as np

# ---------- Utility functions ----------

def euc_dist(p, q):
    """Euclidean distance between two 3D coordinates."""
    return math.sqrt((q[0]-p[0])**2 + (q[1]-p[1])**2 + (q[2]-p[2])**2)

def get_atom_coords(pdb_lines, resname=None, resnum=None, atom_names=None):
    """Extract coordinates filtered by residue name/number and atom names."""
    coords = []
    for line in pdb_lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        res = line[17:20].strip()
        res_id = line[22:26].strip()
        atom = line[12:16].strip()
        if (resname is None or res == resname) and \
           (resnum is None or res_id == str(resnum)) and \
           (atom_names is None or atom in atom_names):
            try:
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append((x, y, z))
            except ValueError:
                pass
    return coords

# ---------- Core function ----------

def add_na_between_carboxylates(pdb_path, res1, res2, offset=0.0):
    """
    Places Na+ at the midpoint between carboxylate oxygens of two residues.
    Optionally offset along the connecting vector by a few Å.
    """
    with open(pdb_path) as f:
        lines = f.readlines()

    res_centers = []
    for resname, resnum in [res1, res2]:
        oxy_names = ["OD1", "OD2"] if resname.startswith("ASP") else ["OE1", "OE2"]
        oxy_coords = get_atom_coords(lines, resname=resname, resnum=resnum, atom_names=oxy_names)
        if not oxy_coords:
            print(f"[WARNING] {resname} {resnum} missing carboxylate atoms.")
            continue
        center = np.mean(np.array(oxy_coords), axis=0)
        res_centers.append(center)
        print(f"[INFO] {resname}{resnum} carboxylate center: {center}")

    if len(res_centers) < 2:
        print("[ERROR] Need two acidic residues with valid oxygens.")
        return lines

    # Midpoint between the two carboxylate centers
    r1, r2 = np.array(res_centers[0]), np.array(res_centers[1])
    midpoint = (r1 + r2) / 2.0

    # Offset along the connecting vector
    vec = r2 - r1
    vec /= np.linalg.norm(vec)
    na_position = midpoint + offset * vec

    cx, cy, cz = na_position
    new_line = (
        f"ATOM   9999  NA   NA I 999    "
        f"{cx:8.3f}{cy:8.3f}{cz:8.3f}  1.00  0.00          NA\n"
    )
    lines.append(new_line)
    print(f"[INFO] Added Na+ between {res1[0]}{res1[1]} and {res2[0]}{res2[1]} "
          f"at ({cx:.3f}, {cy:.3f}, {cz:.3f}) Å.")

    return lines

# ---------- Main ----------

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print(__doc__)
        sys.exit(1)

    pdb_file = sys.argv[1]
    res1 = (sys.argv[2], int(sys.argv[3]))
    res2 = (sys.argv[4], int(sys.argv[5]))
    offset = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0

    new_lines = add_na_between_carboxylates(pdb_file, res1, res2, offset=offset)

    output_file = pdb_file.replace(".pdb", "_with_NA.pdb")
    with open(output_file, "w") as f:
        f.writelines(new_lines)

    print(f"[INFO] Wrote output file: {output_file}")
