"""Step 3: predict phi(z)/trains/loops/tails/adsorbed amount for the real
PSMA-stat-PMMA and PSMA-block-PMMA sequences, given the chi/chis derived in
steps 1-2. Pure sf1d Python API call, no MD -- this part has no blocking
dependency once results/chi.yaml and results/chis.yaml exist.

config["validate"]["sequences"]["stat"/"block"] are still null placeholders
(see pipeline/config.yaml) -- fill in the real composition/degree of
polymerization before this script does anything useful. A stat (random)
sequence is a plain list of "A"/"B" characters in whatever order; sf1d
accepts that directly as an array (see docs/usage.md on the Python API),
not just the block-string CLI syntax.
"""

import sys

import numpy as np
import yaml

sys.path.insert(0, "../../src")  # fall back to the sf1d source tree if not pip-installed
from sf1d import SF1D  # noqa: E402


def main():
    with open(snakemake.input.chi) as fh:
        chi = yaml.safe_load(fh)
    with open(snakemake.input.chis) as fh:
        chis = yaml.safe_load(fh)

    sequences = snakemake.config["validate"]["sequences"]
    results = {}
    for label, seq in sequences.items():
        if seq is None:
            results[label] = {"status": "skipped: no real sequence configured yet"}
            continue
        arr = np.asarray(list(seq))
        m = SF1D(arr, chi=chi, chis=chis, phib=snakemake.config.get("phib", 1e-4),
                  M=snakemake.config.get("M", 80)).solve()
        summary, profiles = m.analyse()
        results[label] = summary

    with open(snakemake.output[0], "w") as fh:
        yaml.safe_dump(results, fh, sort_keys=False)


if __name__ == "__main__":
    main()
