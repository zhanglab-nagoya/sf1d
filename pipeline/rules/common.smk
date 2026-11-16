"""Shared helpers: species table, mdp rendering, GROMACS grompp/mdrun rules.

Species for the bulk (step0/step1) stage:
- HEX:  hexadecane itself (already a complete itp; no polyply gen_params needed)
- ST:   capped stearyl-methacrylate oligomer,  CH3n-(ST)_n-CH3n
- MMA:  capped methyl-methacrylate oligomer,   CH3n-(MMA)_n-CH3n
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


def backbone_bond_path_length(itp_path, n_atoms):
    """BFS shortest path (bond count) from atom 1 to atom n_atoms in an
    itp's [ bonds ] section. Used by vacuum_nvt to size its anti-folding
    pull restraint from the oligomer's own topology instead of a hardcoded
    chain length -- side groups only add tree branches off the backbone,
    never a shortcut, so the shortest path between the two end-cap atoms
    is guaranteed to follow the backbone itself."""
    import collections

    adj = collections.defaultdict(list)
    in_bonds = False
    with open(itp_path) as fh:
        for line in fh:
            s = line.strip()
            if s.startswith("["):
                in_bonds = s.lower().startswith("[ bonds")
                continue
            if not in_bonds or not s or s.startswith(";"):
                continue
            a, b = int(s.split()[0]), int(s.split()[1])
            adj[a].append(b)
            adj[b].append(a)

    dist = {1: 0}
    queue = collections.deque([1])
    while queue:
        u = queue.popleft()
        if u == n_atoms:
            return dist[u]
        for v in adj[u]:
            if v not in dist:
                dist[v] = dist[u] + 1
                queue.append(v)
    raise ValueError(f"{itp_path}: no bonded path from atom 1 to atom {n_atoms}")
