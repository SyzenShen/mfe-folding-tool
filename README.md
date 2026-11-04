# MFE folding tool

A small local web app for practising nucleic acid secondary-structure prediction
(minimum free energy, MFE) for DNA and RNA, e.g. for exercises like
"fold 5 random 40-mers as DNA and RNA and compare free energy and paired bases".

**This is not NUPACK.** It uses [ViennaRNA](https://www.tbi.univie.ac.at/RNA/)
(Turner 2004 parameters for RNA, Mathews 2004 for DNA). Energy parameters and
algorithm details differ from NUPACK, so free energies will not match NUPACK
exactly. The difference has **not been quantified yet**. Use this for practice and
cross-checking, and follow your course instructions for graded values.

![Screenshot: three RNA sequences folded, with table and structure drawings](docs/screenshot.jpg)

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

## Run in Docker / hosting

A `Dockerfile` is included (gunicorn, non-root user, request limits lowered to
20 sequences of up to 300 nt via `MFE_MAX_SEQS` / `MFE_MAX_LEN`). It listens on
`$PORT`, defaulting to 7860.

```bash
docker build -t mfe-tool .
docker run --rm -p 7860:7860 mfe-tool      # then open http://127.0.0.1:7860
```

Built and run locally on an Intel Mac (linux/amd64) and checked against the
non-Docker run. It has **not** been deployed to any hosting service yet.
Notes from the providers' documentation (read 2026-11-04, not tested by me):

- Hugging Face: Docker Spaces need a paid plan to create (PRO for personal
  accounts); free Spaces go to sleep when unused.
- Render: free web services spin down after 15 minutes without traffic and take
  about a minute to wake up.

## Limits

- Single strands only, no multi-strand complexes or concentrations.
- Local default is up to 50 sequences of up to 500 nt per request. Measured on
  an Intel i7 laptop: one 300 nt sequence takes about 0.4 s and one 500 nt
  sequence about 1.5-2 s (RNA, including the drawing). "Both" doubles the work.
- `python app.py` uses Flask's development server on localhost only. For anything
  public use the Docker image (gunicorn). There is no rate limiting or
  authentication.

## Credit and licence

Folding and the partition function are computed with ViennaRNA, Institute for
Theoretical Chemistry, University of Vienna. ViennaRNA's licence allows
research, educational and commercial use without a fee, provided credit is
given; putting it into a commercial product requires contacting its authors.
See the ViennaRNA package for the full text.

The code in this repository is released under the MIT licence (see `LICENSE`).
ViennaRNA is a separate dependency with its own licence.
