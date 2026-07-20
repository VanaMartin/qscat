# CSVE.V00.J00 — N₂ vibrational-excitation cross sections (golden data)

- **Source:** Karel Houfek, time-independent calculation. External to this repo (not from eMoScat).
- **System:** electron–N₂, ²Π_g resonance (LCP model). Initial state v=0, J=0.
- **Format:** 400 rows × 32 whitespace-separated columns, Fortran `E` notation.
  - Column 1: collision energy, **Hartree** (5e-4 … 0.2, step 5e-4).
  - Column 2: elastic / vibrationally-elastic (v=0→0).
  - Columns 3–32: v=0→1, v=0→2, …, v=0→30.
- **Units:** cross sections in **atomic units (bohr²)**.
- Higher-v channels are exactly 0 below their energetic threshold.
- Used as regression anchors for a future time-independent solver (see reference.py).
