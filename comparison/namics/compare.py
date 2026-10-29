"""Cross-check sf1d against Namics (github.com/leermakers/Namics).

Each case is run in both codes with the same lattice, sequence, chi values and
Silberberg adsorption energies chis.  In Namics a wall is a frozen segment
type 'Si' (and 'Si2' for the second wall) just outside the system, so
chi_{x,Si} = -chis_x / lambda1 = -6 chis_x (solvent: chi_{W,Si} = 0).

Usage:  python compare.py /path/to/namics/bin/namics
Writes results.csv and gap_comparison.png next to this script.
"""

import pathlib
import re
import subprocess
import sys

import numpy as np

from sf1d import SF1D, parse_sequence

LAM1 = 1.0 / 6.0
HERE = pathlib.Path(__file__).resolve().parent


def namics_input(seq, chi, chis, phib, M, walls=1, chis2=None, theta=None):
    """Namics input text.  chi keys like 'AS' (S = solvent), chis like {'A': 2.0}."""
    chis2 = chis if chis2 is None else chis2
    types = sorted(set(parse_sequence(seq)))
    comp = "".join(f"({t}){n}" for t, n in re.findall(r"([A-RT-Za-z])(\d+)", seq))
    name = lambda t: "W" if t == "S" else t
    lines = [
        f"lat : flat : n_layers : {M}",
        "lat : flat : lattice_type : simple_cubic",
        "lat : flat : geometry : planar",
        "lat : flat : gradients : 1",
        "lat : flat : lowerbound : surface",
        f"lat : flat : upperbound : {'surface' if walls == 2 else 'mirror'}",
        "mon : Si : freedom : frozen",
        "mon : Si : frozen_range : lowerbound",
        "mon : W : freedom : free",
    ]
    if walls == 2:
        lines += ["mon : Si2 : freedom : frozen", "mon : Si2 : frozen_range : upperbound"]
    for t in types:
        lines.append(f"mon : {t} : freedom : free")
        if chis.get(t, 0.0):
            lines.append(f"mon : {t} : chi_Si : {-chis[t] / LAM1:.12g}")
        if walls == 2 and chis2.get(t, 0.0):
            lines.append(f"mon : {t} : chi_Si2 : {-chis2[t] / LAM1:.12g}")
    for key, val in chi.items():
        lines.append(f"mon : {name(key[0])} : chi_{name(key[1])} : {val:.12g}")
    lines += ["mol : water : composition : (W)1", "mol : water : freedom : solvent",
              f"mol : pol : composition : {comp}"]
    lines += ["mol : pol : freedom : free", f"mol : pol : phibulk : {phib:.12g}",
              "newton : isaac : method : pseudohessian",
              "newton : isaac : iterationlimit : 5000",
              "newton : isaac : tolerance : 1e-11",
              "output : pro : append : false", "output : kal : append : false",
              "kal : sys : NN : grand_potential", "kal : sys : NN : free_energy",
              "kal : mol : pol : theta", "kal : mol : pol : phibulk",
              "kal : mol : pol : mu", "kal : mol : pol : n",
              "kal : mol : water : mu", "kal : mol : water : n"]
    lines += [f"pro : mon : {t} : phi" for t in types] + ["pro : mon : W : phi", "start"]
    if theta is not None:
        # restricted: start from the full-equilibrium solution at the same gap
        # (Namics' default guess fails for compressed restricted layers);
        # the second run overwrites the output files
        lines += ["sys : NN : initial_guess : previous_result",
                  "mol : pol : freedom : restricted", f"mol : pol : theta : {theta:.15g}", "start"]
    lines.append("")
    return "\n".join(lines)


def run_namics(binary, tag, **case):
    (HERE / "inputs").mkdir(exist_ok=True)
    (HERE / "output").mkdir(exist_ok=True)
    inp = HERE / "inputs" / f"{tag}.in"
    inp.write_text(namics_input(**case))
    log = subprocess.run([binary, f"inputs/{tag}.in"], cwd=HERE, capture_output=True, text=True).stdout
    (HERE / "output" / f"{tag}.log").write_text(log)
    if log.count("Problem solved") != log.count("Problem nr") or "negative phibulk" in log:
        raise RuntimeError(f"Namics did not converge for {tag}; see output/{tag}.log")
    # several 'start' blocks number the profile files: <tag>_1.pro, <tag>_2.pro, ...
    pros = sorted((HERE / "output").glob(f"{tag}_[0-9]*.pro")) or [HERE / "output" / f"{tag}.pro"]
    pro = np.genfromtxt(pros[-1], names=True)
    kal = np.atleast_1d(np.genfromtxt(HERE / "output" / f"{tag}.kal", names=True))  # one row per run
    return pro, kal


def compare(binary, tag, seq, chi, chis, phib, M, walls=1, chis2=None, theta=None):
    """Run one case in both codes and compare profiles, amounts and free energies.

    Full equilibrium: Namics' grand potential vs sf1d's Omega.
    Restricted: Namics' Helmholtz F is referenced to its own bulk, which it
    rescales with the effective polymer bulk fraction.  Referenced to the true
    reservoir (as in sf1d), the free energy is
        A = F - n_S mu_S(reservoir) - n_p mu_p(reservoir),
    with the reservoir chemical potentials taken from the full-equilibrium
    run that precedes the restricted run in the same Namics input."""
    pro, kal = run_namics(binary, tag, seq=seq, chi=chi, chis=chis, phib=phib, M=M,
                          walls=walls, chis2=chis2, theta=theta)
    types = sorted(set(parse_sequence(seq)))
    m = SF1D(seq, chi=chi, chis=chis, chis2=chis2, phib=phib, M=M, walls=walls, theta=theta).solve()
    phi_sf = m.phi[[m.idx[t] for t in types]]
    phi_nm = np.array([np.atleast_1d(pro[f"mon_{t}_phi"]) for t in types])
    last = kal[-1]
    if theta is None:
        F_nm = last["sys_NN_grand_potential"]
    else:
        F_nm = (last["sys_NN_free_energy"] - last["mol_water_n"] * kal[0]["mol_water_mu"]
                - last["mol_pol_n"] * kal[0]["mol_pol_mu"])
    return dict(case=tag, max_dphi=float(np.max(np.abs(phi_sf - phi_nm))),
                theta_sf1d=float(phi_sf.sum()), theta_namics=float(last["mol_pol_theta"]),
                F_sf1d=m.free_energy(), F_namics=float(F_nm))


AB = dict(seq="A10B90", chi={"AS": 0.8, "BS": 0.4}, chis={"A": 2.0}, phib=1e-4)
CASES = [
    dict(tag="homo_athermal", seq="A100", chi={}, chis={"A": 1.0}, phib=1e-4, M=80),
    dict(tag="homo_theta_solvent", seq="A100", chi={"AS": 0.5}, chis={"A": 1.0}, phib=1e-4, M=80),
    dict(tag="anchor_buoy", M=80, **AB),
    dict(tag="anchor_buoy_slit_D20", M=20, walls=2, **AB),
    dict(tag="anchor_buoy_slit_D30_inert_wall2", M=30, walls=2, chis2={}, **AB),
]


def all_cases():
    """CASES plus a compressed restricted slit: the amount at D=40 squeezed into D=10."""
    theta = SF1D(**AB, M=40, walls=2).solve().analyse()[0]["theta_total"]
    return CASES + [dict(tag="anchor_buoy_slit_D10_restricted", M=10, walls=2, theta=theta, **AB)]


def gap_comparison(binary, Ds=(6, 8, 10, 12, 15, 20, 25, 30, 35, 40, 50, 60)):
    """dF(D) = F(D) - 2 gamma in full and restricted equilibrium from both codes
    (anchor-buoy A10B90; restricted amount frozen at its full value at max(Ds))."""
    Ds = sorted(Ds, reverse=True)
    g = compare(binary, "gap_single_wall", M=80, **AB)
    rows, theta_fix = [], None
    for mode in ("full", "restricted"):
        for D in Ds:
            try:
                r = compare(binary, f"gap_{mode}_D{D}", M=D, walls=2,
                            theta=None if mode == "full" else theta_fix, **AB)
            except RuntimeError as e:
                print(f"  {mode} D={D}: {e}")
                continue
            if mode == "full" and D == Ds[0]:
                theta_fix = r["theta_namics"]
            rows.append((mode, D, r["F_sf1d"] - 2 * g["F_sf1d"], r["F_namics"] - 2 * g["F_namics"], r["max_dphi"]))
    return rows


if __name__ == "__main__":
    binary = sys.argv[1] if len(sys.argv) > 1 else "namics"
    out = ["case,max_dphi,theta_sf1d,theta_namics,F_sf1d,F_namics"]
    print(f"{'case':34s} {'max|dphi|':>10s} {'theta sf1d':>12s} {'theta Namics':>12s} {'F sf1d':>12s} {'F Namics':>12s}")
    for case in all_cases():
        r = compare(binary, **case)
        print(f"{r['case']:34s} {r['max_dphi']:10.2e} {r['theta_sf1d']:12.6f} {r['theta_namics']:12.6f} "
              f"{r['F_sf1d']:12.6f} {r['F_namics']:12.6f}")
        out.append(",".join(str(r[k]) for k in ("case", "max_dphi", "theta_sf1d", "theta_namics", "F_sf1d", "F_namics")))
    (HERE / "results.csv").write_text("\n".join(out) + "\n")

    rows = gap_comparison(binary)
    lines = ["mode,D,dF_sf1d,dF_namics,max_dphi"] + [",".join(map(str, r)) for r in rows]
    (HERE / "gap_comparison.csv").write_text("\n".join(lines) + "\n")
    print("gap scan: max |dF_sf1d - dF_namics| =",
          f"{max(abs(r[2] - r[3]) for r in rows):.2e}, max |dphi| = {max(r[4] for r in rows):.2e}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(5.5, 4))
    for mode, c in (("full", "C0"), ("restricted", "C1")):
        d = [r for r in rows if r[0] == mode]
        ax.plot([r[1] for r in d], [r[2] for r in d], "-", color=c, label=f"sf1d, {mode}")
        ax.plot([r[1] for r in d], [r[3] for r in d], "o", mfc="none", color=c, label=f"Namics, {mode}")
    ax.axhline(0, color="grey", lw=0.8)
    ax.set(xlabel="gap D (layers)", ylabel="dF (kT per site)", title="A10B90 between two walls: sf1d vs Namics")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(HERE / "gap_comparison.png", dpi=150)
