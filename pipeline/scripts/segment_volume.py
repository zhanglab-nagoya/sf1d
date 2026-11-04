"""Step 0: lattice segment volume a^3 (from the hexadecane NPT density) and
how many lattice segments one ST/MMA repeat unit maps to (from the
oligomer-melt NPT densities and molar masses read straight out of the
force field files already in forcefield/ -- nothing hardcoded).

Run as a snakemake `script:` -- `snakemake` is injected (input/output/params).
"""

import re

import numpy as np
import yaml

N_A = 6.02214076e23  # 1/mol


def read_xvg_density(path):
    """Average density (kg/m^3) from a gmx energy -o density.xvg file."""
    vals = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "@")):
                continue
            vals.append(float(line.split()[1]))
    return float(np.mean(vals))


def itp_molar_mass(ff_path):
    """Sum of atom masses (g/mol) in a GROMACS [ atoms ] block."""
    mass = 0.0
    in_atoms = False
    with open(ff_path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("["):
                in_atoms = s.lower().startswith("[ atoms")
                continue
            if not in_atoms or not s or s.startswith(";"):
                continue
            fields = s.split()
            # nr type resnr residue atom cgnr charge mass
            if len(fields) >= 8:
                mass += float(fields[7])
    return mass


def main():
    ff_dir = snakemake.config["forcefield_dir"]
    n_oligo = snakemake.params.n_oligo

    density = {
        sp: read_xvg_density(f"results/bulk/{sp}/density.xvg")
        for sp in snakemake.params.species
    }

    # V_molecule[nm^3] = (M[g/mol] / N_A[1/mol]) / rho[kg/m^3] * 1e24
    m_hex = itp_molar_mass(f"{ff_dir}/solvent/hexadecane.itp")
    a3_nm3 = (m_hex / N_A) / density["HEX"] * 1e24

    result = {
        "a3_nm3": a3_nm3,
        "hexadecane_molar_mass_g_per_mol": m_hex,
        "hexadecane_density_kg_per_m3": density["HEX"],
        "segments_per_repeat_unit": {},
    }
    for sp in ("ST", "MMA"):
        m_monomer = itp_molar_mass(f"{ff_dir}/monomers/PAMA_{sp}.oplsaa.LigParGen.ff")
        # Partial molar volume of one repeat unit, approximated as the whole
        # oligomer's mass density (caps are a small fraction of a 10-mer's
        # mass; this approximation gets worse for shorter oligomer_length).
        v_repeat_nm3 = (m_monomer / N_A) / density[sp] * 1e24
        result["segments_per_repeat_unit"][sp] = v_repeat_nm3 / a3_nm3
        result[f"{sp}_oligomer_density_kg_per_m3"] = density[sp]
        result[f"{sp}_monomer_molar_mass_g_per_mol"] = m_monomer

    with open(snakemake.output[0], "w") as fh:
        yaml.safe_dump(result, fh, sort_keys=False)


if __name__ == "__main__":
    main()
