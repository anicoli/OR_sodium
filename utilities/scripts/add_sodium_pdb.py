#!/usr/bin/env python3
"""
add_na_between_acids.py
=======================
Place a sodium ion (Na+) at the midpoint between the carboxylate groups
of two acidic residues (ASP or GLU) in a PDB file.

Usage
-----
    python add_na_between_acids.py input.pdb RES1 NUM1 RES2 NUM2 [OFFSET]

Arguments
---------
input.pdb  : input PDB file
RES1 NUM1  : residue name and number of the first acidic residue  (e.g. ASP 70)
RES2 NUM2  : residue name and number of the second acidic residue (e.g. ASP 111)
OFFSET     : optional displacement along the inter-residue vector in Å (default: 0.0)

Example
-------
    python add_na_between_acids.py protein.pdb ASP 70 ASP 111 0.5

Output
------
    input_with_NA.pdb  –  copy of the input PDB with the Na+ ATOM record appended.
    If the output file already exists, the script will exit without overwriting.
"""

import os
import sys
import numpy as np


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

# Carboxylate oxygen atom names for each supported residue type
CARBOXYLATE_OXYGENS = {
    "ASP": ["OD1", "OD2"],
    "GLU": ["OE1", "OE2"],
}


def get_carboxylate_center(pdb_lines: list, resname: str, resnum: int) -> np.ndarray:
    """
    Return the geometric centre of the carboxylate oxygens for a given residue.

    Parameters
    ----------
    pdb_lines : list of str  – lines read from the PDB file
    resname   : str          – residue name (ASP or GLU)
    resnum    : int          – residue sequence number

    Returns
    -------
    np.ndarray of shape (3,)  – centre coordinates in Å

    Raises
    ------
    ValueError  if the residue type is not supported or no oxygens are found
    """
    resname = resname.upper()
    if resname not in CARBOXYLATE_OXYGENS:
        raise ValueError(
            f"Unsupported residue '{resname}'. Supported types: "
            + ", ".join(CARBOXYLATE_OXYGENS)
        )

    target_atoms = set(CARBOXYLATE_OXYGENS[resname])
    coords = []

    for line in pdb_lines:
        if not line.startswith(("ATOM", "HETATM")):
            continue
        if line[17:20].strip() != resname:
            continue
        if line[22:26].strip() != str(resnum):
            continue
        if line[12:16].strip() not in target_atoms:
            continue
        try:
            coords.append([
                float(line[30:38]),
                float(line[38:46]),
                float(line[46:54]),
            ])
        except ValueError:
            pass  # skip malformed coordinate fields

    if not coords:
        raise ValueError(
            f"No carboxylate oxygens ({', '.join(target_atoms)}) found "
            f"for {resname} {resnum}."
        )

    center = np.mean(np.array(coords), axis=0)
    print(f"[INFO] {resname} {resnum} carboxylate centre: "
          f"({center[0]:.3f}, {center[1]:.3f}, {center[2]:.3f}) Å  "
          f"({len(coords)} oxygens used)")
    return center


def next_serial(pdb_lines: list) -> int:
    """
    Return the next available ATOM serial number (max existing + 1, capped at 99999).

    Parameters
    ----------
    pdb_lines : list of str

    Returns
    -------
    int
    """
    max_serial = 0
    for line in pdb_lines:
        if line.startswith(("ATOM", "HETATM")):
            try:
                max_serial = max(max_serial, int(line[6:11]))
            except ValueError:
                pass
    return min(max_serial + 1, 99999)


def format_na_record(serial: int, x: float, y: float, z: float) -> str:
    """
    Format a PDB ATOM record for a sodium ion following the PDB fixed-width standard.

    Columns (1-based):
      1-6   record type
      7-11  serial
      13-16 atom name
      18-20 residue name
      22    chain ID
      23-26 residue sequence number
      31-38 x  (8.3f)
      39-46 y  (8.3f)
      47-54 z  (8.3f)
      55-60 occupancy
      61-66 B-factor
      77-78 element symbol

    Parameters
    ----------
    serial : int
    x, y, z : float  – coordinates in Å

    Returns
    -------
    str  – properly formatted PDB ATOM record (no trailing newline)
    """
    return (
        f"{'ATOM':<6}{serial:5d}  {'NA':<4} {'NA':>3} {'I':1}{999:4d}    "
        f"{x:8.3f}{y:8.3f}{z:8.3f}{1.00:6.2f}{0.00:6.2f}          {'NA':>2}"
    )


# ---------------------------------------------------------------------------
# Core function
# ---------------------------------------------------------------------------

def add_na_between_carboxylates(
    pdb_path: str,
    res1: tuple,
    res2: tuple,
    offset: float = 0.0,
) -> list:
    """
    Add a Na+ ion between the carboxylate groups of two acidic residues.

    The ion is placed at the midpoint of the two carboxylate centres, with an
    optional displacement along the inter-residue axis.

    Parameters
    ----------
    pdb_path : str             – path to the input PDB file
    res1     : (str, int)      – (resname, resnum) of the first residue
    res2     : (str, int)      – (resname, resnum) of the second residue
    offset   : float           – displacement along the inter-residue vector (Å)

    Returns
    -------
    list of str  – PDB lines including the new Na+ ATOM record
    """
    with open(pdb_path) as fh:
        lines = fh.readlines()

    # Locate carboxylate centres
    c1 = get_carboxylate_center(lines, res1[0], res1[1])
    c2 = get_carboxylate_center(lines, res2[0], res2[1])

    # Place Na+ at the midpoint, optionally displaced along the connecting vector
    midpoint = (c1 + c2) / 2.0
    direction = c2 - c1
    norm = np.linalg.norm(direction)
    if norm < 1e-6:
        raise ValueError("The two carboxylate centres are at the same position.")
    direction /= norm

    na_pos = midpoint + offset * direction
    cx, cy, cz = na_pos

    # Report inter-centre distance as a sanity check (ideal Na–O ~2.4 Å)
    dist = np.linalg.norm(c2 - c1)
    print(f"[INFO] Distance between carboxylate centres: {dist:.2f} Å")
    if dist > 10.0:
        print("[WARNING] Centres are more than 10 Å apart — verify residue selection.")

    serial  = next_serial(lines)
    na_line = format_na_record(serial, cx, cy, cz) + "\n"
    lines.append(na_line)

    print(f"[INFO] Na+ added (serial {serial}) at "
          f"({cx:.3f}, {cy:.3f}, {cz:.3f}) Å "
          f"between {res1[0]} {res1[1]} and {res2[0]} {res2[1]}.")

    return lines


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if len(sys.argv) < 6:
        print(__doc__)
        sys.exit(1)

    pdb_file = sys.argv[1]
    res1     = (sys.argv[2].upper(), int(sys.argv[3]))
    res2     = (sys.argv[4].upper(), int(sys.argv[5]))
    offset   = float(sys.argv[6]) if len(sys.argv) > 6 else 0.0

    # Validate input file
    if not os.path.isfile(pdb_file):
        print(f"[ERROR] Input PDB not found: {pdb_file}", file=sys.stderr)
        sys.exit(1)

    # Prevent silent overwrite
    output_file = pdb_file.replace(".pdb", "_with_NA.pdb")
    if os.path.isfile(output_file):
        print(f"[ERROR] Output file already exists: {output_file}\n"
              f"        Remove or rename it before running.", file=sys.stderr)
        sys.exit(1)

    try:
        new_lines = add_na_between_carboxylates(pdb_file, res1, res2, offset=offset)
    except ValueError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)

    with open(output_file, "w") as fh:
        fh.writelines(new_lines)

    print(f"[INFO] Output written to: {output_file}")
