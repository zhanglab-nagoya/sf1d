# Theory and implementation

`sf1d` implements the Scheutjens–Fleer (SF) lattice self-consistent field theory for
linear (co)polymers in a monomeric solvent, near one flat wall or between two parallel
walls. Concentrations vary only in the direction normal to the walls (z); within a layer
the mean-field (Bragg–Williams) approximation is used.

## Lattice

Simple cubic lattice, layers `z = 1..M`. A site in layer z has a fraction
λ₀ = 4/6 of its neighbours in the same layer and λ₁ = 1/6 in each adjacent layer.
The site average of a quantity f is

    <f>(z) = λ₁ f(z−1) + λ₀ f(z) + λ₁ f(z+1)

All species (every polymer segment and the solvent) occupy one site of volume a³.

**Boundaries.** z = 0 is always a wall. With one wall, z = M+1 is a mirror
(f(M+1) = f(M)), so M must be large enough for the profile to reach the bulk.
With two walls, z = M+1 is a second wall and the gap width is D = M layers.

## Segment potentials

For segment type x (polymer types and solvent S), in units of kT:

    u_x(z) = u'(z) + Σ_y χ_xy (<φ_y>(z) − φ_y^b) − χs_x δ(z,1) − χs2_x δ(z,M)

- u'(z) is the Lagrange field that enforces incompressibility, Σ_x φ_x(z) = 1.
- χ_xy are Flory–Huggins parameters (symmetric, zero diagonal).
- χs_x is the adsorption exchange parameter of segment x relative to solvent
  (Silberberg convention): χs > 0 means x is preferred over solvent in the wall layer.
  The solvent has χs = 0 by definition. χs2 applies to the second wall (default: χs).
- In the mixing term, a wall's neighbour slot is counted as bulk (φ_y = φ_y^b), so a
  wall exerts no force beyond χs. With this convention the critical adsorption energy
  of long chains in an athermal solvent is χs_c = −ln(1 − λ₁) = ln(6/5) ≈ 0.182, the
  textbook cubic-lattice value.

The bulk has u = 0 for every species.

## Chain statistics

Free-segment weights are G_x(z) = exp(−u_x(z)). For a chain with segment sequence
t(1..N), the forward and backward propagators are

    G(z, 1|1) = G_t(1)(z),   G(z, s|1) = G_t(s)(z) <G(·, s−1|1)>(z)
    G(z, N|N) = G_t(N)(z),   G(z, s|N) = G_t(s)(z) <G(·, s+1|N)>(z)

with G = 0 at walls. The volume fraction of segment s is

    φ(z, s) = C G(z, s|1) G(z, s|N) / G_t(s)(z)

and φ_x(z) sums φ(z, s) over all s of type x. The solvent is monomeric:
φ_S(z) = φ_S^b G_S(z).

**Normalisation.**
- *Full equilibrium* (default): the gap exchanges chains with a bulk solution of
  polymer volume fraction φ_b, so C = φ_b / N.
- *Restricted equilibrium* (`theta`): the total polymer amount θ = Σ_z φ_p(z) in the
  gap is fixed, so C = θ / (N Q) with Q = Σ_z G(z, N|1) the single-chain partition
  function. The solvent still exchanges with the bulk.

## Solver

The unknowns are u_x(z) for all types and layers. The residual is

    r_x(z) = u'_x(z) − mean_y u'_y(z) + ln Σ_y φ_y(z),     u'_x = u_x − (mixing + surface terms)

which vanishes when all types see the same Lagrange field and the layer is full.
Using ln Σφ (rather than 1 − 1/Σφ) keeps Newton steps bounded far from the solution.

Picard and Anderson mixing are unstable here for long adsorbing chains, because φ
depends on u roughly as exp(N·u). `SF1D.solve` uses damped Newton instead: a
forward-difference Jacobian, evaluated as one batched set of perturbations, and a
backtracking line search. It typically converges in 10–30 steps from u = 0.
`gap_scan` warm-starts each gap from the previous one and retries cold if that fails.

## Free energy

The grand potential per lattice site of wall area, in kT and relative to bulk, is

    Ω = Σ_z [ −(φ_p − φ_b)/N − (φ_S − φ_S^b) − Σ_x φ_x u_x + ½ Σ_xy χ_xy δφ_x <δφ_y>₀ ]
        − Σ_x χs_x φ_x(1) − Σ_x χs2_x φ_x(M)            (second term: two walls only)

with δφ = φ − φ^b and <·>₀ the site average with wall slots set to zero. Its stationary
point reproduces the SCF equations above; the tests check this by finite differences
(∂Ω/∂χs_A = −φ_A in the wall layers; ∂Ω/∂χ_AS equals the explicit derivative of the
mixing term).

- One wall: Ω is the interfacial tension γ (kT/a²).
- Restricted equilibrium: the free energy is A = Ω + (θ/N) ln(C N / φ_b), which equals
  Ω when θ takes its full-equilibrium value.
- Two walls: the interaction free energy is ΔF(D) = Ω(D) − γ₁ − γ₂, and the normal
  pressure is Π(D) = −dΔF/dD (central differences over the scanned gaps), in kT/a³.
  Π > 0 is repulsive.

In `gap_scan`, restricted equilibrium freezes θ at its full-equilibrium value at the
largest gap, mimicking compression that is fast compared with chain exchange. The scan
stops once θ would exceed 90 % of the gap.

## Adsorbed-chain analysis

A chain is *adsorbed* if at least one segment is in a wall layer. Propagators computed
with the wall layers forbidden (G = 0 there) give the free chains; the difference gives
the adsorbed ones. Following the standard SF decomposition:

- **trains**: adsorbed segments in a wall layer,
- **tails**: adsorbed segments whose path to one chain end touches no wall,
- **loops**: the remaining adsorbed segments (between two wall contacts, including bridges).

With two walls, chains touching wall 1, wall 2, or both are counted by inclusion–exclusion
over propagators that forbid layer 1, layer M, or both. Chains touching both walls are
**bridges**.

## Limitations

- **Mean field.** Excluded-volume correlations are averaged away. In good solvents the
  scaling exponents are wrong (e.g. a dilute adsorbed layer decays as z⁻² instead of
  the de Gennes z⁻⁴/³). Trends, adsorbed amounts and qualitative structure are usually
  right.
- **Equilibrium.** The theory assumes the layer has equilibrated. Strongly adsorbing
  segments (several kT) can freeze in non-equilibrium structures; see
  O'Shaughnessy & Vavylonis (2003).
- **Coarse chemistry.** One site size for every species, one χs per segment type, a
  homogeneous wall. Side-chain crowding of comb polymers, chain stiffness, specific
  bonding and surface heterogeneity are not represented.
- **No bulk microphases.** A 1D calculation cannot represent micelles in the bulk.
  For long, poorly soluble blocks (high χ_AS), check that the bulk is a single phase.
- **Statics only.** No dynamics, flow or friction.

## References

- J. M. H. M. Scheutjens and G. J. Fleer, *J. Phys. Chem.* **83**, 1619 (1979); **84**, 178 (1980).
- J. M. H. M. Scheutjens and G. J. Fleer, Interaction between two adsorbed polymer layers,
  *Macromolecules* **18**, 1882 (1985).
- O. A. Evers, J. M. H. M. Scheutjens and G. J. Fleer, Statistical thermodynamics of block
  copolymer adsorption. 1, *Macromolecules* **23**, 5221 (1990).
- G. J. Fleer, M. A. Cohen Stuart, J. M. H. M. Scheutjens, T. Cosgrove and B. Vincent,
  *Polymers at Interfaces* (Chapman and Hall, 1993).
- P. G. de Gennes, *Macromolecules* **14**, 1637 (1981).
- B. O'Shaughnessy and D. Vavylonis, *Phys. Rev. Lett.* **90**, 056103 (2003);
  *Eur. Phys. J. E* **11**, 213 (2003).
