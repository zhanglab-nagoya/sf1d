"""Step 3: validate the derived chi/chis together against full-chain MD.

sf1d_predict (below) is real and runnable once results/chi.yaml and
results/chis.yaml exist (steps 1-2) and config["validate"]["sequences"] is
filled in with the actual PSMA-stat-PMMA / PSMA-block-PMMA composition.

The MD side of this comparison (full real-sequence chain + hexadecane +
Fe2O3 substrate, production run, measured adsorbed amount/profile) is design
only, same status as step2_chis.smk -- it's the same build/convert/merge
machinery as step 2, just with a real multi-segment chain instead of a
single fragment, so it should reuse whatever step 2 ends up being rather
than inventing a second version. Not scaffolded here; revisit once step 2
is working.
"""


rule sf1d_predict:
    input:
        chi="results/chi.yaml",
        chis="results/chis.yaml",
    output:
        "results/sf1d_predictions.yaml",
    conda:
        "../envs/analysis.yaml"
    script:
        "../scripts/sf1d_predict.py"
