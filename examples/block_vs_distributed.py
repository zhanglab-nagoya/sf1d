"""Same anchor content (10 % A), different placement along the chain:
one A block vs. ten single A segments spread out."""

import json
import sys

from sf1d import SF1D

params = dict(chis={"A": 2.0}, chi={"AS": 0.8, "BS": 0.4}, phib=1e-4, M=80)
cases = {"block A10B90": "A10B90", "distributed (A1B9)x10": "A1B9" * 10}
keys = ("theta_adsorbed", "adsorbed_chains_per_site", "fraction_trains", "fraction_loops",
        "fraction_tails", "rms_distance_to_nearest_wall_adsorbed")

results = {}
for name, seq in cases.items():
    summary, _ = SF1D(seq, **params).solve().analyse()
    results[name] = {k: round(summary[k], 4) for k in keys}
    print(f"{name:24s} " + "  ".join(f"{k}={results[name][k]}" for k in keys))

out = sys.argv[1] if len(sys.argv) > 1 else "block_vs_distributed.json"
with open(out, "w") as fh:
    json.dump({"parameters": params, "results": results}, fh, indent=2)
