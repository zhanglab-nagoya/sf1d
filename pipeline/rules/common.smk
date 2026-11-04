"""Shared helpers: species table, mdp rendering, GROMACS grompp/mdrun rules.

Species for the bulk (step0/step1) stage:
- HEX:  hexadecane itself (already a complete itp; no polyply gen_params needed)
- ST:   capped stearyl-methacrylate oligomer,  CH3-(ST)_n-CH3n
- MMA:  capped methyl-methacrylate oligomer,   CH3-(MMA)_n-CH3n
"""

FF_DIR = config["forcefield_dir"]
N_OLIGO = config["bulk"]["oligomer_length"]

BULK_SPECIES = ["HEX", "ST", "MMA"]


def molname(species):
    return "HEX" if species == "HEX" else f"{species}{N_OLIGO}"


def render_mdp(template_path, out_path, **fmt):
    with open(template_path) as fh:
        text = fh.read()
    with open(out_path, "w") as fh:
        fh.write(text.format(**fmt))
