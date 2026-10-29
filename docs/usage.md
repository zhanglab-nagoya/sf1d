# Usage

## Parameters

| Option | Meaning |
|---|---|
| `--seq` | Block sequence of segment types, e.g. `A10B90` or `A1B9A1B9`. Any letter except `S` (the solvent). |
| `--chi XY=v` | Flory–Huggins parameter between types X and Y, e.g. `AS=0.8`. Repeatable; unset pairs are 0. |
| `--chis X=v` | Adsorption exchange parameter of type X at wall 1, relative to solvent. Repeatable; unset types are 0. |
| `--chis2 X=v` | Same for wall 2 (two walls only). Default: same as `--chis`. Use e.g. `--chis2 A=0` for a non-adsorbing counter-surface. |
| `--phib` | Bulk polymer volume fraction (default 1e-4). |
| `--M` | Number of lattice layers. One wall: domain size (must reach the bulk). Two walls: gap width D. |
| `--walls {1,2}` | One wall with a mirror at z = M+1, or two walls. |
| `--theta` | Restricted equilibrium: fixed total polymer amount in the gap (segments per site). |
| `--scan DMIN DMAX` | Two-wall gap scan from DMAX down to DMIN (implies two walls). |
| `--step` | Gap step for `--scan` (default 1). |
| `--mode {full,restricted,both}` | Equilibrium mode for `--scan` (default both). |
| `--out` | Output file prefix. |
| `--plot` | Also write a PNG (needs matplotlib). |

All energies are in kT, lengths in lattice units a, amounts in segments per lattice site
of wall area.

## Single profile

```bash
sf1d --seq A10B90 --chis A=2.0 --chi AS=0.8 --chi BS=0.4 --phib 1e-4 --M 80 --out run --plot
sf1d --seq A10B90 --chis A=2.0 --chi AS=0.8 --chi BS=0.4 --walls 2 --M 20 --out slit
```

Writes:

- `run_summary.json` (also printed):
  `free_energy` (Ω, or A in restricted equilibrium), `theta_total`, `theta_excess`,
  `theta_adsorbed`, `fraction_trains/loops/tails` of adsorbed segments,
  `adsorbed_chains_per_site`, `rms_distance_to_nearest_wall_adsorbed`,
  per-type `theta_adsorbed_X` and `train_X`, and for two walls
  `bridging_chains_per_site` and `bridging_fraction_of_adsorbed_chains`.
- `run_profiles.csv`: columns `z, total_polymer, adsorbed, free, trains, loops, tails, phi_<type>...`.
- `run_profiles.png` with `--plot`.

## Gap scan

```bash
sf1d --seq A10B90 --chis A=2.0 --chi AS=0.8 --chi BS=0.4 --phib 1e-4 \
     --scan 4 60 --mode both --out gap --plot
```

Writes `gap_scan_full.csv` and/or `gap_scan_restricted.csv` with columns

| Column | Meaning |
|---|---|
| `D` | Gap width (layers) |
| `dF` | Interaction free energy ΔF = Ω(D) − γ₁ − γ₂ (kT per site) |
| `Pi` | Normal pressure −dΔF/dD (kT/a³), > 0 repulsive |
| `theta` | Total polymer amount in the gap |
| `bridges` | Bridging chains per site |
| `bridge_frac` | Bridging fraction of adsorbed chains |

and `gap_scan.png` with `--plot`. Choose DMAX large enough that ΔF(DMAX) ≈ 0 (the
scan warns otherwise); the restricted amount is taken from the full-equilibrium
state at DMAX.

**Units.** kT/a³ is 4.1 MPa for a = 1 nm and 33 MPa for a = 0.5 nm (T ≈ 298 K).

## Python API

```python
from sf1d import SF1D, gap_scan

m = SF1D("A10B90", chis={"A": 2.0}, chi={"AS": 0.8, "BS": 0.4}, phib=1e-4, M=80).solve()
summary, profiles = m.analyse()       # dicts; profiles hold numpy arrays
m.phi                                 # (n_types, M) volume fractions, order m.types
m.free_energy()                       # gamma (one wall), Omega or A (two walls)

slit = SF1D("A10B90", chis={"A": 2.0}, chis2={"A": 0.0}, chi={"AS": 0.8, "BS": 0.4},
            phib=1e-4, M=20, walls=2).solve()

res = gap_scan("A10B90", chi={"AS": 0.8, "BS": 0.4}, chis={"A": 2.0}, chis2=None,
               phib=1e-4, Ds=range(4, 61), mode="both")
res["restricted"]["D"], res["restricted"]["Pi"]
```

`solve(u0=...)` accepts a previous `m.u` as a warm start (same M).

## Using the same parameters in Namics

`chis` is the Silberberg χs of the standard Scheutjens–Fleer model, so values carry
over to Namics and the literature directly: `--chis A=v` corresponds to a wall segment
`Si` with `mon : A : chi_Si : −6v`. The full mapping, and how to compare
restricted-equilibrium free energies, is in
[comparison/namics/README.md](../comparison/namics/README.md).

## Practical notes

- Choose the segment size first, e.g. one methacrylate monomer with its side chain, or
  one Kuhn segment. Every species, including the solvent, occupies one site of that size.
- Treat χ and χs as effective parameters. They can be fitted to MD, for example χs from
  monomer adsorption PMFs and χ from solvation free energies or single-chain R_g.
- Very long chains (N ≳ 1000) at the tightest gaps can defeat Newton. The scan then
  retries from a cold start and otherwise stops with a message.
