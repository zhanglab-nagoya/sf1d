"""Command-line interface: single profiles (one or two walls) and gap scans."""

import argparse
import json

import numpy as np

from . import core
from .core import SF1D, gap_scan


def kv_pairs(items):
    out = {}
    for item in items or []:
        k, v = item.split("=")
        out[k] = float(v)
    return out


def main():
    p = argparse.ArgumentParser(prog="sf1d", description=core.__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--seq", default="A10B90", help="block sequence, e.g. A10B90")
    p.add_argument("--chi", action="append", help="Flory-Huggins pair, e.g. AS=0.8 (repeatable)")
    p.add_argument("--chis", action="append", help="surface exchange, e.g. A=2.0 (repeatable)")
    p.add_argument("--phib", type=float, default=1e-4, help="bulk polymer volume fraction")
    p.add_argument("--M", type=int, default=80, help="number of lattice layers (= gap D for two walls)")
    p.add_argument("--walls", type=int, default=1, choices=(1, 2), help="one wall + mirror, or two walls")
    p.add_argument("--chis2", action="append", help="second-wall surface exchange (default: same as --chis)")
    p.add_argument("--theta", type=float, help="restricted equilibrium: fixed polymer amount in the gap")
    p.add_argument("--scan", type=int, nargs=2, metavar=("DMIN", "DMAX"), help="two-wall gap scan")
    p.add_argument("--step", type=int, default=1, help="gap step for --scan")
    p.add_argument("--mode", default="both", choices=("full", "restricted", "both"), help="equilibrium for --scan")
    p.add_argument("--out", default="sf1d_run", help="output prefix")
    p.add_argument("--plot", action="store_true", help="write a PNG of the profiles")
    a = p.parse_args()

    chis2 = kv_pairs(a.chis2) if a.chis2 else None
    if a.scan:
        run_scan(a, chis2)
        return

    model = SF1D(a.seq, chi=kv_pairs(a.chi), chis=kv_pairs(a.chis), phib=a.phib, M=a.M,
                 walls=a.walls, chis2=chis2, theta=a.theta).solve()
    summary, prof = model.analyse()
    summary.update({"seq": a.seq, "chi": kv_pairs(a.chi), "chis": kv_pairs(a.chis), "M": a.M})
    if a.walls == 2:
        summary["chis2"] = chis2 if chis2 is not None else kv_pairs(a.chis)

    with open(f"{a.out}_summary.json", "w") as fh:
        json.dump(summary, fh, indent=2)
    keys = list(prof)
    np.savetxt(f"{a.out}_profiles.csv", np.column_stack([prof[k] for k in keys]),
               delimiter=",", header=",".join(keys), comments="")
    print(json.dumps(summary, indent=2))

    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
        z = prof["z"]
        for t in model.poly_types:
            ax1.semilogy(z, prof[f"phi_{t}"], label=f"segment {t}")
        ax1.semilogy(z, prof["total_polymer"], "k--", lw=1, label="all polymer")
        ax1.axhline(a.phib, color="grey", lw=0.8, ls=":", label="bulk")
        ax1.set(xlabel="layer z", ylabel="volume fraction", title="Profiles by segment type")
        ax1.legend(frameon=False)
        inner = slice(1, None) if a.walls == 1 else slice(1, -1)
        for k in ("loops", "tails", "free"):
            ax2.semilogy(z[inner], prof[k][inner], label=k)
        wl = [0] if a.walls == 1 else [0, -1]
        ax2.semilogy(z[wl], prof["trains"][wl], "ko", label="trains (wall layers)")
        ax2.set(xlabel="layer z", ylabel="volume fraction", title="Adsorbed-chain decomposition")
        ax2.set_ylim(bottom=max(1e-8, a.phib * 1e-3))
        ax2.legend(frameon=False)
        fig.suptitle(f"{a.seq}, chis={kv_pairs(a.chis)}, chi={kv_pairs(a.chi)}, phib={a.phib:g}", fontsize=9)
        fig.tight_layout()
        fig.savefig(f"{a.out}_profiles.png", dpi=150)


def run_scan(a, chis2):
    Ds = np.arange(a.scan[0], a.scan[1] + 1, a.step)
    res = gap_scan(a.seq, kv_pairs(a.chi), kv_pairs(a.chis), chis2, a.phib, Ds, mode=a.mode)
    modes = [m for m in ("full", "restricted") if m in res]
    print(f"gamma_1 = {res['gamma'][0]:.6f}, gamma_2 = {res['gamma'][1]:.6f} kT/a^2; "
          f"restricted theta = {res['theta_fixed']:.4f}")
    for md in modes:
        r = res[md]
        cols = ["D", "dF", "Pi", "theta", "bridges", "bridge_frac"]
        np.savetxt(f"{a.out}_scan_{md}.csv", np.column_stack([r[c] for c in cols]),
                   delimiter=",", header=",".join(cols), comments="")
        i = np.argmin(r["dF"])
        print(f"{md:10s}: min dF = {r['dF'][i]:.4f} kT/a^2 at D = {r['D'][i]:.0f}; "
              f"Pi at D = {r['D'][-1]:.0f}: {r['Pi'][-1]:.4f} kT/a^3")
    if a.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, axs = plt.subplots(1, 3, figsize=(13, 4))
        for md in modes:
            r = res[md]
            axs[0].plot(r["D"], r["dF"], "o-", ms=3, label=md)
            axs[1].plot(r["D"], r["Pi"], "o-", ms=3, label=md)
            axs[2].plot(r["D"], r["bridge_frac"], "o-", ms=3, label=md)
        for ax in axs[:2]:
            ax.axhline(0, color="grey", lw=0.8)
        axs[0].set(xlabel="gap D (layers)", ylabel="dF (kT per site)", title="Interaction free energy")
        axs[1].set(xlabel="gap D (layers)", ylabel="Pi (kT / a^3)", title="Normal pressure (>0 repulsive)")
        axs[1].set_yscale("symlog", linthresh=1e-4)
        axs[2].set(xlabel="gap D (layers)", ylabel="fraction", title="Bridging fraction of adsorbed chains")
        for ax in axs:
            ax.legend(frameon=False)
        fig.suptitle(f"{a.seq}, chis={kv_pairs(a.chis)}, chis2={chis2 or 'same'}, chi={kv_pairs(a.chi)}, "
                     f"phib={a.phib:g}", fontsize=9)
        fig.tight_layout()
        fig.savefig(f"{a.out}_scan.png", dpi=150)


if __name__ == "__main__":
    main()
