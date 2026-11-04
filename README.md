# MFE folding tool

A small local web app for practising nucleic acid secondary-structure prediction
(minimum free energy, MFE) for DNA and RNA, e.g. for exercises like
"fold 5 random 40-mers as DNA and RNA and compare free energy and paired bases".

**This is not NUPACK.** It uses [ViennaRNA](https://www.tbi.univie.ac.at/RNA/)
(Turner 2004 parameters for RNA, Mathews 2004 for DNA). Energy parameters and
algorithm details differ from NUPACK, so free energies will not match NUPACK
exactly. The difference has **not been quantified yet**. Use this for practice and
cross-checking, and follow your course instructions for graded values.

## Features

- Several sequences at once (one per line, optionally `name SEQUENCE`); DNA, RNA or both.
- Temperature setting (default 37 °C).
- MFE structure (dot-bracket), ΔG in kcal/mol, number of base pairs, number of
  paired bases (2 x pairs), and averages over all sequences.
- Structure drawings with three layouts (Puzzler, Circular, NAView), bases
  coloured by their probability of being paired in the thermodynamic ensemble,
  scroll to zoom, drag to pan, double-click to reset.

## Run

Needs Python 3.8+.

```bash
python3 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python app.py
```

Then open <http://127.0.0.1:5050>. (Port 5050 is used because macOS AirPlay
Receiver often occupies 5000.)

## Test

```bash
.venv/bin/python test_app.py
```

The tests check the code runs consistently (parsing, T/U handling, layouts,
probability ranges). They do not validate the energy values against NUPACK.

## Limits

- Single strands only, no multi-strand complexes or concentrations.
- Up to 50 sequences of up to 500 nt per request; long sequences were not benchmarked.
- Uses Flask's development server and binds to localhost only. It is not set up
  for public hosting.

## Credit

Folding and the partition function are computed with ViennaRNA, Institute for
Theoretical Chemistry, University of Vienna. ViennaRNA's licence allows
research and educational use without a fee, provided credit is given. See the
ViennaRNA package for the full text.
