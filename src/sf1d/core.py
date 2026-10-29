#!/usr/bin/env python3
"""
sf1d.core -- minimal 1D Scheutjens-Fleer self-consistent field (SCF) solver
for a linear copolymer adsorbing from solution onto one flat wall, or into
the slit between two parallel walls.

Lattice and conventions
-----------------------
* Simple cubic lattice, layers z = 1..M parallel to the wall.
  A segment in layer z has a fraction lambda0 = 4/6 of its neighbours in z
  and lambda1 = 1/6 in each of z-1 and z+1.
* z = 0 is an impenetrable wall.  With one wall (default), z = M+1 is a
  mirror (reflecting) boundary, so M must be large enough that the profile
  reaches the bulk.  With two walls (--walls 2), z = M+1 is a second wall
  with its own adsorption parameters (--chis2, default = --chis) and the gap
  width is D = M lattice layers.
* Components: polymer segment types (any letters, e.g. 'A', 'B') plus the
  monomeric solvent 'S'.  The chain sequence is given as a block string,
  e.g. "A10B90" = 10 anchor segments followed by 90 buoy segments.
* Full equilibrium (default): the polymer is in equilibrium with a bulk
  solution of volume fraction phib.  Potentials are measured relative to the
  bulk (u = 0 in the bulk), so the normalisation is C = phib / N.
  Restricted equilibrium (--theta): the total amount of polymer between the
  walls, theta (segments per site), is fixed; the solvent still exchanges
  with the bulk.  The normalisation then follows from theta.
* Segment potential (units of kT):
      u_x(z) = u'(z) + sum_y chi_xy (<phi_y>(z) - phi_y^b)
               - chis_x delta(z,1) - chis2_x delta(z,M)   [second term: 2 walls]
  u'(z) is the Lagrange field enforcing incompressibility, sum_x phi_x = 1.
  chis_x is the Silberberg adsorption exchange parameter of segment x
  relative to solvent (solvent chis = 0; chis > 0 means adsorbing).
  For the mixing term, a wall's neighbour slot is counted as bulk, so a
  wall exerts no force beyond chis.  With this convention the
  critical adsorption energy for long chains is chis_c = -ln(1 - lambda1)
  = ln(6/5) ~= 0.182, the textbook cubic-lattice value.

* Grand potential per lattice site of wall area (kT units), relative to bulk:
      Omega = sum_z [ -sum_i (phi_i - phi_i^b)/N_i - sum_x phi_x u_x
                      + 1/2 sum_xy chi_xy dphi_x <dphi_y>_0 ]
              - sum_x chis_x phi_x(1) [- sum_x chis2_x phi_x(M)]
  with dphi = phi - phi^b and <.>_0 the site average with wall slots = 0.
  Its stationary point reproduces the SCF equations above.  For restricted
  equilibrium the free energy is A = Omega + (theta/N) ln(C N / phib).
  With one wall Omega is the interfacial tension gamma; with two walls the
  interaction free energy is dF(D) = Omega(D) - gamma_1 - gamma_2, and the
  normal pressure is Pi = -d dF / dD (kT per a^3).

Outputs
-------
Volume-fraction profiles per segment type, the adsorbed-chain profile split
into trains / loops / tails (standard SF decomposition), and a summary with
adsorbed amounts, bridging chains (two walls) and the free energy.  --scan
computes dF(D) and Pi(D) over a range of gaps, in full and/or restricted
equilibrium.

Usage
-----
    sf1d --seq A10B90 --chis A=2.0 --chi AS=0.8 --chi BS=0.4 \
                   --phib 1e-4 --M 80 --out run1
    sf1d --seq A10B90 --chis A=2.0 --chi AS=0.8 --chi BS=0.4 \
                   --walls 2 --scan 4 60 --mode both --out gap --plot
"""

import re

import numpy as np

LAM0, LAM1 = 4.0 / 6.0, 1.0 / 6.0


# --------------------------------------------------------------------------- #
# Model setup
# --------------------------------------------------------------------------- #
def parse_sequence(seq):
    """'A10B90' -> np.array(['A']*10 + ['B']*90)."""
    blocks = re.findall(r"([A-RT-Za-z])(\d+)", seq)  # 'S' reserved for solvent
    if not blocks or "".join(f"{t}{n}" for t, n in blocks) != seq:
        raise ValueError(f"Cannot parse sequence '{seq}' (use e.g. A10B90; 'S' is the solvent).")
    return np.array([t for t, n in blocks for _ in range(int(n))])


class SF1D:
    def __init__(self, seq, chi=None, chis=None, phib=1e-4, M=80, walls=1, chis2=None, theta=None):
        self.chain = parse_sequence(seq) if isinstance(seq, str) else np.asarray(seq)
        self.N = len(self.chain)
        self.poly_types = sorted(set(self.chain))
        self.types = self.poly_types + ["S"]
        self.idx = {t: i for i, t in enumerate(self.types)}
        self.M = M
        self.phib = phib
        if walls not in (1, 2):
            raise ValueError("walls must be 1 or 2")
        self.walls = walls
        self.theta = theta  # None: full equilibrium; float: restricted equilibrium
        n = len(self.types)

        # Flory-Huggins matrix (symmetric, zero diagonal)
        self.chi = np.zeros((n, n))
        for key, val in (chi or {}).items():
            a, b = key[0], key[1]
            i, j = self.idx[a], self.idx[b]
            self.chi[i, j] = self.chi[j, i] = val

        # Surface exchange parameters (solvent = 0 by definition)
        self.chis = np.zeros(n)
        for t, val in (chis or {}).items():
            self.chis[self.idx[t]] = val
        self.chis2 = np.zeros(n)
        for t, val in (chis if chis2 is None else chis2 or {}).items():
            self.chis2[self.idx[t]] = val

        # Bulk composition
        self.phibulk = np.zeros(n)
        for t in self.poly_types:
            self.phibulk[self.idx[t]] = phib * np.mean(self.chain == t)
        self.phibulk[self.idx["S"]] = 1.0 - phib

        self.seg_type = np.array([self.idx[t] for t in self.chain])

    # ----------------------------------------------------------------------- #
    # All array functions accept optional leading batch dimensions (...),
    # so the Newton Jacobian can be evaluated for many perturbations at once.
    def avg(self, f, wall_value):
        """Site average <f>(z); wall slots take wall_value, mirror at M+1 (1 wall)."""
        pad = np.broadcast_to(wall_value, f.shape[:-1] + (1,))
        fm = np.concatenate((pad, f[..., :-1]), axis=-1)
        fp = np.concatenate((f[..., 1:], pad if self.walls == 2 else f[..., -1:]), axis=-1)
        return LAM1 * fm + LAM0 * f + LAM1 * fp

    def step(self, G):
        """Propagator step <G>(z): G = 0 at walls, mirror at M+1 (1 wall)."""
        return self.avg(G, 0.0)

    def propagators(self, Gw, mask=None):
        """Forward/backward chain propagators.  Gw: (..., ntypes, M) free-segment
        weights.  mask (M,) multiplies every step, e.g. to forbid layer 1."""
        Gseg = Gw[..., self.seg_type, :]  # (..., N, M)
        if mask is not None:
            Gseg = Gseg * mask
        Gf = np.empty_like(Gseg)
        Gb = np.empty_like(Gseg)
        N = self.N
        Gf[..., 0, :] = Gseg[..., 0, :]
        for s in range(1, N):
            Gf[..., s, :] = Gseg[..., s, :] * self.step(Gf[..., s - 1, :])
        Gb[..., N - 1, :] = Gseg[..., N - 1, :]
        for s in range(N - 2, -1, -1):
            Gb[..., s, :] = Gseg[..., s, :] * self.step(Gb[..., s + 1, :])
        return Gf, Gb, Gseg

    def densities(self, u):
        """Volume fractions (..., ntypes, M) from potentials u (..., ntypes, M)."""
        Gw = np.exp(-u)
        Gf, Gb, Gseg = self.propagators(Gw)
        phi_seg = self.norm(Gf) * Gf * Gb / Gseg  # (..., N, M)
        phi = np.zeros_like(u)
        for t in self.poly_types:
            k = self.idx[t]
            phi[..., k, :] = phi_seg[..., self.seg_type == k, :].sum(axis=-2)
        kS = self.idx["S"]
        phi[..., kS, :] = self.phibulk[kS] * Gw[..., kS, :]
        return phi

    def norm(self, Gf):
        """Polymer normalisation C, shape (..., 1, 1): phib/N in full equilibrium,
        theta/(N Q) in restricted equilibrium (Q = single-chain partition function)."""
        if self.theta is None:
            return np.full(Gf.shape[:-2] + (1, 1), self.phib / self.N)
        Q = Gf[..., -1, :].sum(axis=-1)
        return (self.theta / (self.N * Q))[..., None, None]

    def interactions(self, phi):
        """Mixing + surface part of u_x(z) (everything except u')."""
        avgphi = self.avg(phi, self.phibulk[:, None])
        inter = np.einsum("xy,...ym->...xm", self.chi, avgphi - self.phibulk[:, None])
        inter[..., 0] -= self.chis
        if self.walls == 2:
            inter[..., -1] -= self.chis2
        return inter

    def residual(self, x):
        """x: (..., ntypes*M) potentials -> residual of the same shape."""
        n = len(self.types)
        u = x.reshape(x.shape[:-1] + (n, self.M))
        phi = self.densities(u)
        uprime = u - self.interactions(phi)
        # incompressibility as ln(sum phi): keeps Newton steps bounded far from the solution
        res = uprime - uprime.mean(axis=-2, keepdims=True) + np.log(phi.sum(axis=-2, keepdims=True))
        return res.reshape(x.shape)

    def jacobian(self, x, r, h=1e-7, batch=64):
        """Forward-difference Jacobian, evaluated in batches of perturbations."""
        n = x.size
        J = np.empty((n, n))
        for i0 in range(0, n, batch):
            idx = np.arange(i0, min(i0 + batch, n))
            X = np.repeat(x[None, :], idx.size, axis=0)
            X[np.arange(idx.size), idx] += h
            J[:, idx] = ((self.residual(X) - r) / h).T
        return J

    # ----------------------------------------------------------------------- #
    def solve(self, u0=None, tol=1e-10, maxiter=100):
        """Damped Newton (finite-difference Jacobian, backtracking line search).
        Picard/Anderson mixing is unstable here for long adsorbing chains,
        because phi depends on u roughly as exp(N * u)."""
        n = len(self.types)
        x = np.zeros(n * self.M) if u0 is None else np.asarray(u0, float).ravel().copy()
        r = self.residual(x)
        for it in range(maxiter):
            if np.max(np.abs(r)) < tol:
                break
            dx = np.linalg.lstsq(self.jacobian(x, r), -r, rcond=None)[0]
            lam, f0 = 1.0, r @ r
            while lam > 1e-6:
                with np.errstate(over="ignore", invalid="ignore", divide="ignore"):
                    rn = self.residual(x + lam * dx)
                    accept = np.all(np.isfinite(rn)) and rn @ rn < (1 - 1e-4 * lam) * f0
                if accept:
                    break
                lam *= 0.5
            else:
                raise RuntimeError(f"line search failed at iteration {it} (residual {np.max(np.abs(r)):.2e})")
            x, r = x + lam * dx, rn
        else:
            raise RuntimeError(f"SCF did not converge in {maxiter} Newton steps (residual {np.max(np.abs(r)):.2e})")
        self.iterations = it
        self.u = x.reshape(n, self.M)
        self.phi = self.densities(self.u)
        self.C = float(self.norm(self.propagators(np.exp(-self.u))[0]).squeeze())
        return self

    def free_energy(self):
        """Grand potential Omega (full equilibrium) or free energy A (restricted),
        per lattice site of wall area, in kT, relative to the bulk."""
        kS = self.idx["S"]
        phi_p = self.phi[[self.idx[t] for t in self.poly_types]].sum(axis=0)
        ideal = -(phi_p - self.phib) / self.N - (self.phi[kS] - self.phibulk[kS])
        dphi = self.phi - self.phibulk[:, None]
        mix = 0.5 * np.einsum("xy,xm,ym->m", self.chi, dphi, self.avg(dphi, 0.0))
        omega = (ideal - (self.phi * self.u).sum(axis=0) + mix).sum()
        omega -= self.chis @ self.phi[:, 0]
        if self.walls == 2:
            omega -= self.chis2 @ self.phi[:, -1]
        if self.theta is not None:
            omega += self.theta / self.N * np.log(self.C * self.N / self.phib)
        return float(omega)

    # ----------------------------------------------------------------------- #
    def analyse(self):
        """Adsorbed/free split, train/loop/tail decomposition and, for two
        walls, the number of bridging chains (touching both walls)."""
        Gw = np.exp(-self.u)
        C = self.C
        Gf, Gb, Gseg = self.propagators(Gw)
        wall_layers = [0] if self.walls == 1 else [0, self.M - 1]
        mask = np.ones(self.M)
        mask[wall_layers] = 0.0  # free chains never touch a wall layer
        Ff, Fb, _ = self.propagators(Gw, mask=mask)
        Af, Ab = Gf - Ff, Gb - Fb  # weights of walks that have touched layer 1

        phi_tot = C * Gf * Gb / Gseg
        with np.errstate(invalid="ignore", divide="ignore"):
            phi_free = np.where(Gseg * mask > 0, C * Ff * Fb / Gseg, 0.0)
            phi_tail = np.where(Gseg * mask > 0, C * (Ff * Ab + Af * Fb) / Gseg, 0.0)
        phi_ads = phi_tot - phi_free

        tot = phi_tot.sum(axis=0)
        ads = phi_ads.sum(axis=0)
        tail = phi_tail.sum(axis=0)
        train = np.zeros(self.M)
        train[wall_layers] = ads[wall_layers]
        loop = ads - train - tail
        loop[wall_layers] = 0.0

        z = np.arange(1, self.M + 1)
        zwall = z - 0.5 if self.walls == 1 else np.minimum(z - 0.5, self.M + 0.5 - z)
        theta_ads = ads.sum()
        summary = {
            "walls": self.walls,
            "mode": "full" if self.theta is None else "restricted",
            "sequence_length_N": int(self.N),
            "phib": self.phib,
            "free_energy": self.free_energy(),
            "theta_total": float(tot.sum()),
            "theta_excess": float((tot - self.phib).sum()),
            "theta_adsorbed": float(theta_ads),
            "fraction_trains": float(train.sum() / theta_ads) if theta_ads > 0 else 0.0,
            "fraction_loops": float(loop.sum() / theta_ads) if theta_ads > 0 else 0.0,
            "fraction_tails": float(tail.sum() / theta_ads) if theta_ads > 0 else 0.0,
            "adsorbed_chains_per_site": float(theta_ads / self.N),
            "rms_distance_to_nearest_wall_adsorbed": float(np.sqrt((ads * zwall ** 2).sum() / theta_ads))
            if theta_ads > 0 else 0.0,
            "max_incompressibility_error": float(np.max(np.abs(self.phi.sum(axis=0) - 1))),
        }
        # adsorbed amount per segment type, and its bound fraction
        for t in self.poly_types:
            k = self.idx[t]
            ads_t = phi_ads[self.seg_type == k].sum(axis=0)
            summary[f"theta_adsorbed_{t}"] = float(ads_t.sum())
            summary[f"train_{t}"] = float(ads_t[wall_layers].sum())
        if self.walls == 2:
            # chains per site touching wall 1 / wall 2 / both, by inclusion-exclusion
            Q = lambda m: self.propagators(Gw, mask=m)[0][-1].sum()
            m1, m2 = np.ones(self.M), np.ones(self.M)
            m1[0], m2[-1] = 0.0, 0.0
            Qall, Qnot1, Qnot2, Qfree = Gf[-1].sum(), Q(m1), Q(m2), Ff[-1].sum()
            bridges = C * (Qall - Qnot1 - Qnot2 + Qfree)
            summary["bridging_chains_per_site"] = float(bridges)
            summary["bridging_fraction_of_adsorbed_chains"] = float(bridges / (C * (Qall - Qfree)))
        profiles = {"z": z, "total_polymer": tot, "adsorbed": ads, "free": tot - ads,
                    "trains": train, "loops": loop, "tails": tail}
        for t in self.types:
            profiles[f"phi_{t}"] = self.phi[self.idx[t]]
        return summary, profiles


# --------------------------------------------------------------------------- #
# Two-wall gap scan
# --------------------------------------------------------------------------- #
def remap(u, M_new):
    """Warm start for a new gap: keep the wall regions, drop/duplicate the middle."""
    M = u.shape[-1]
    k = (M_new + 1) // 2
    right = u[:, M - (M_new - k):] if M_new - k > 0 else u[:, :0]
    if M_new > M:  # widen: pad the middle with the central value
        mid = np.repeat(u[:, M // 2:M // 2 + 1], M_new - M, axis=1)
        return np.concatenate((u[:, :M // 2], mid, u[:, M // 2:]), axis=1)
    return np.concatenate((u[:, :k], right), axis=1)


def gap_scan(seq, chi, chis, chis2, phib, Ds, mode="both", verbose=True):
    """Interaction free energy dF(D) and pressure Pi(D) between two walls.
    Full equilibrium: chains exchange with the bulk.  Restricted: the amount
    of polymer in the gap is frozen at its full-equilibrium value at max(Ds).
    Returns {mode: dict of arrays}."""
    Ds = np.sort(np.asarray(Ds, int))[::-1]  # compress from wide to narrow
    common = dict(chi=chi, phib=phib)
    Msingle = max(80, int(Ds[0]))
    gamma1 = SF1D(seq, chis=chis, M=Msingle, **common).solve().free_energy()
    gamma2 = gamma1 if chis2 is None else SF1D(seq, chis=chis2, M=Msingle, **common).solve().free_energy()
    modes = ["full", "restricted"] if mode == "both" else [mode]
    out, theta_fix = {}, None
    for md in modes:
        rows, u = [], None
        for D in Ds:
            theta = None if md == "full" else theta_fix
            if theta is not None and theta > 0.9 * D:
                if verbose:
                    print(f"  {md}: stop at D={D}, fixed theta={theta:.2f} would exceed 90% of the gap")
                break
            m = SF1D(seq, chis=chis, chis2=chis2, M=int(D), walls=2, theta=theta, **common)
            try:
                m.solve(u0=None if u is None else remap(u, int(D)))
            except RuntimeError:
                try:  # warm start can land outside Newton's basin: retry cold
                    m.solve(maxiter=300)
                except RuntimeError as e:
                    if verbose:
                        print(f"  {md}: no convergence at D={D} ({e}); stopping")
                    break
            u = m.u
            s, _ = m.analyse()
            if md == "full" and D == Ds[0]:
                theta_fix = s["theta_total"]
            rows.append((D, s["free_energy"] - gamma1 - gamma2, s["theta_total"],
                         s["bridging_chains_per_site"], s["bridging_fraction_of_adsorbed_chains"]))
        r = np.array(rows)
        res = dict(D=r[:, 0], dF=r[:, 1], theta=r[:, 2], bridges=r[:, 3], bridge_frac=r[:, 4])
        res["Pi"] = -np.gradient(res["dF"], res["D"]) if len(r) > 1 else np.zeros(len(r))
        out[md] = res
        if verbose and abs(res["dF"][0]) > 1e-3 * abs(gamma1 + gamma2):
            print(f"  warning ({md}): dF(Dmax) = {res['dF'][0]:.2e} is not ~0; increase Dmax")
    out["gamma"] = (gamma1, gamma2)
    out["theta_fixed"] = theta_fix
    return out


# --------------------------------------------------------------------------- #
# Analytic reference
# --------------------------------------------------------------------------- #
def critical_chis_transfer_matrix(M=400):
    """Largest eigenvalue of the single-chain transfer matrix crosses 1 at
    chis_c; returns the numerically located value (dilute, athermal)."""
    from scipy.optimize import brentq

    def lam_max(chis):
        T = np.diag(np.full(M, LAM0)) + np.diag(np.full(M - 1, LAM1), 1) + np.diag(np.full(M - 1, LAM1), -1)
        T[-1, -1] += LAM1  # mirror boundary
        T[0] *= np.exp(chis)
        return np.max(np.linalg.eigvals(T).real)

    return brentq(lambda c: lam_max(c) - 1.0, 0.05, 0.5)
