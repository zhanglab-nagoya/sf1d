#!/usr/bin/env bash
# Regenerate every example output in examples/output/.  Requires `pip install -e ".[plot]"`.
set -euo pipefail
cd "$(dirname "$0")"
mkdir -p output
AB="--chis A=2.0 --chi AS=0.8 --chi BS=0.4 --phib 1e-4"

echo "1) Anchor-buoy A10B90 on a single wall"
sf1d --seq A10B90 $AB --M 80 --out output/anchor_buoy --plot > /dev/null

echo "2) Block vs. distributed anchors"
python block_vs_distributed.py output/block_vs_distributed.json

echo "3) Anchor-buoy A10B90 between two walls: gap scan"
sf1d --seq A10B90 $AB --scan 4 60 --mode both --out output/gap --plot

echo "4) Validation: N=1000 homopolymer gap scan (Scheutjens-Fleer 1985 trends)"
sf1d --seq A1000 --chis A=1.0 --chi AS=0.5 --phib 1e-4 --scan 3 40 --mode both --out output/homopolymer_N1000 --plot
