#!/usr/bin/env python3
import sys
import re
from decimal import Decimal, getcontext

getcontext().prec = 28  # plenty for charge arithmetic

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} input.itp output.itp")
    sys.exit(1)

inp, outp = sys.argv[1], sys.argv[2]

PRECISION = 8
EPS = Decimal(f"1e-{PRECISION}")

atom_re = re.compile(
    r"""
    ^(\s*\d+\s+\S+\s+\d+\s+\S+\s+\S+\s+\d+\s+)
    ([+-]?\d+\.\d+)
    (\s+.*)
    """,
    re.VERBOSE,
)

# ---------- First pass: read charges ----------
charges = []
lines = []
in_atoms = False

with open(inp) as f:
    for line in f:
        stripped = line.strip()

        if stripped.startswith("[ atoms ]"):
            in_atoms = True
            lines.append((line, None))
            continue

        if in_atoms:
            if stripped.startswith("["):
                in_atoms = False
                lines.append((line, None))
                continue

            if stripped and not stripped.startswith(";"):
                m = atom_re.match(line.rstrip("\n"))
                if m:
                    q = Decimal(m.group(2))
                    charges.append(q)
                    lines.append((line, q))
                else:
                    lines.append((line, None))
            else:
                lines.append((line, None))
        else:
            lines.append((line, None))

N = len(charges)
Q = sum(charges)

delta = Q / Decimal(N)

print(f"Original total charge: {Q:+.10e}")
print(f"Atoms in repeat unit: {N}")
print(f"Primary correction per atom: {-delta:+.10e}")

# ---------- Primary correction + quantization ----------
corrected = []
for q in charges:
    q_new = q - delta
    q_q = q_new.quantize(EPS)
    corrected.append(q_q)

Q2 = sum(corrected)
residual = -Q2  # what we still need to add

steps = int((residual / EPS).to_integral_value())

print(f"Residual after quantization: {Q2:+.10e}")
print(f"Distributing {steps} × {EPS:+.1e}")

# ---------- Distribute residual ----------
sign = 1 if steps > 0 else -1
steps = abs(steps)

for i in range(steps):
    corrected[i] += sign * EPS

assert sum(corrected) == Decimal("0"), "Charge renormalization failed!"

# ---------- Write output ----------
atom_idx = 0
in_atoms = False

with open(outp, "w") as fout:
    for line, q in lines:
        if q is None:
            fout.write(line)
            continue

        stripped = line.strip()
        if stripped.startswith("[ atoms ]"):
            in_atoms = True
            fout.write(line)
            continue

        if in_atoms and stripped.startswith("["):
            in_atoms = False
            fout.write(line)
            continue

        m = atom_re.match(line.rstrip("\n"))
        if m:
            q_fmt = f"{corrected[atom_idx]:.{PRECISION}f}"
            fout.write(m.group(1) + q_fmt + m.group(3) + "\n")
            atom_idx += 1
        else:
            fout.write(line)

print("Done. Charges are exactly neutral at written precision.")
