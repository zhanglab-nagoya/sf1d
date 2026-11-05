"""Step 1: cohesive energy density -> Hildebrand solubility parameters ->
chi_AS, chi_BS (segment_x-solvent Flory-Huggins parameters, A=ST, B=MMA,
S=hexadecane). See forcefield/chi_chis_from_md.md step 1 ("primary, cheap"
route) for the derivation this implements:

    CED_x = (U_vac,x - U_liq,x) / V_molar,x   [per-mole cohesive energy
                                                density of species x's own
                                                melt]
    delta_x = sqrt(CED_x)                     [Hildebrand solubility param]
    chi_xS = V_ref * (delta_x - delta_S)^2 / (R*T)   [V_ref = a^3, the
                                                       solvent-sized lattice
                                                       segment volume from
                                                       step 0]

U_vac,x is the single-molecule vacuum-reference potential energy (rules
vacuum_em/vacuum_nvt); U_liq,x is the per-molecule potential energy in the
melt (rule bulk_potential, same NPT production window as bulk_density).
CED is an intensive quantity -- using the whole ST10/MMA10 oligomer as the
"molecule" (rather than per repeat-unit) gives the same delta_x as long as
U and V are both computed per-oligomer, which is what happens here.

Run as a snakemake `script:` -- `snakemake` is injected (input/output/params).
"""

import math

import numpy as np
import yaml

N_A = 6.02214076e23  # 1/mol
R = 8.314462618  # J/(mol K)


def read_xvg_avg(path):
    """Average of a gmx energy -o *.xvg file's second column."""
    vals = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(("#", "@")):
                continue
            vals.append(float(line.split()[1]))
    return float(np.mean(vals))


def itp_molar_mass(itp_path):
    """Sum of atom masses (g/mol) in a GROMACS [ atoms ] block."""
    mass = 0.0
    in_atoms = False
    with open(itp_path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("["):
                in_atoms = s.lower().startswith("[ atoms")
                continue
            if not in_atoms or not s or s.startswith(";"):
                continue
            fields = s.split()
            if len(fields) >= 8:
                mass += float(fields[7])
    return mass


def parse_n_chains(top_path):
    """Actual molecule count from a topol.top's [ molecules ] section (not
    the requested config count -- see bulk_top's docstring on why packer
    output can differ)."""
    with open(top_path) as fh:
        lines = fh.readlines()
    for i, line in enumerate(lines):
        if line.strip().lower() == "[ molecules ]":
            for data_line in lines[i + 1 :]:
                s = data_line.strip()
                if not s or s.startswith(";"):
                    continue
                return int(s.split()[1])
    raise ValueError(f"no [ molecules ] section found in {top_path}")


def main():
    species = snakemake.params.species  # ["HEX", "ST", "MMA"]
    T = snakemake.config["temperature_K"]

    with open(snakemake.input.segment_volume) as fh:
        seg_vol = yaml.safe_load(fh)

    density_key = {
        "HEX": "hexadecane_density_kg_per_m3",
        "ST": "ST_oligomer_density_kg_per_m3",
        "MMA": "MMA_oligomer_density_kg_per_m3",
    }

    delta = {}  # Hildebrand solubility parameter, Pa^0.5
    ced = {}  # cohesive energy density, J/m^3
    details = {}
    for sp, top_path, itp_path, liq_path, vac_path in zip(
        species,
        snakemake.input.topol,
        snakemake.input.itp,
        snakemake.input.liquid_potential,
        snakemake.input.vacuum_potential,
    ):
        n_chains = parse_n_chains(top_path)
        M = itp_molar_mass(itp_path)  # g/mol, whole oligomer/HEX molecule
        density = seg_vol[density_key[sp]]  # kg/m^3

        v_molar = (M * 1e-3) / density  # m^3/mol

        u_liq_total = read_xvg_avg(liq_path)  # kJ/mol, whole box
        u_liq_per_mol = u_liq_total / n_chains  # kJ/mol, per molecule
        u_vac = read_xvg_avg(vac_path)  # kJ/mol, n=1

        delta_e_coh = (u_vac - u_liq_per_mol) * 1000.0  # J/mol
        if delta_e_coh <= 0:
            raise ValueError(
                f"{sp}: vacuum potential ({u_vac:.1f} kJ/mol) is not above "
                f"the liquid's per-molecule potential ({u_liq_per_mol:.1f} "
                "kJ/mol) -- cohesive energy density would be negative, "
                "which means the melt isn't actually more stable than an "
                "isolated molecule. Check the NPT/vacuum runs before "
                "trusting this number."
            )

        ced[sp] = delta_e_coh / v_molar  # J/m^3
        delta[sp] = math.sqrt(ced[sp])  # Pa^0.5
        details[sp] = {
            "n_chains": n_chains,
            "molar_mass_g_per_mol": M,
            "density_kg_per_m3": density,
            "liquid_potential_kJ_per_mol": u_liq_per_mol,
            "vacuum_potential_kJ_per_mol": u_vac,
            "ced_J_per_m3": ced[sp],
            "delta_MPa0.5": delta[sp] / 1000.0,
        }

    a3_nm3 = seg_vol["a3_nm3"]
    v_ref = a3_nm3 * 1e-27 * N_A  # m^3/mol, the solvent-sized lattice segment

    chi_as = v_ref * (delta["ST"] - delta["HEX"]) ** 2 / (R * T)
    chi_bs = v_ref * (delta["MMA"] - delta["HEX"]) ** 2 / (R * T)

    result = {
        "chi_AS": chi_as,
        "chi_BS": chi_bs,
        "temperature_K": T,
        "v_ref_m3_per_mol": v_ref,
        "species": details,
    }

    with open(snakemake.output[0], "w") as fh:
        yaml.safe_dump(result, fh, sort_keys=False)


if __name__ == "__main__":
    main()
