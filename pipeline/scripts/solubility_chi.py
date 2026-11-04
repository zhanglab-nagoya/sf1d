"""Step 1: cohesive energy density -> Hildebrand solubility parameters ->
chi_AS, chi_BS.

NOT IMPLEMENTED. This needs the vacuum single-molecule reference energy
(rule `vacuum_ref`, also not implemented -- see its docstring in
rules/step0_step1_bulk.smk) to split each melt's U_total into U_intra +
U_inter:

    CED_x = -U_inter_x / V_x
    delta_x = sqrt(CED_x)
    chi_xy = a^3 * (delta_x - delta_y)^2 / (R * T)

Fails loudly on purpose rather than writing fabricated numbers to
results/chi.yaml.
"""

raise NotImplementedError(
    "solubility_chi: waiting on rule `vacuum_ref` (see its docstring). "
    "Do not stub this out with approximate numbers -- implement vacuum_ref "
    "first, then come back here."
)
