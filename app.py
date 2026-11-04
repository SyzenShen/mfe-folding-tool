"""Minimal MFE folding tool for teaching exercises.

Engine: ViennaRNA (https://www.tbi.univie.ac.at/RNA/), Institute for Theoretical
Chemistry, University of Vienna. Free energies can differ from NUPACK's.

Run:  .venv/bin/python app.py     ->  http://127.0.0.1:5050
"""
import math
import os
import re
import threading

import RNA
from flask import Flask, render_template_string, request
from markupsafe import Markup

app = Flask(__name__)

# ViennaRNA loads energy parameters and plot settings globally, so serialize calls.
_lock = threading.Lock()

# Request limits; hosted deployments lower these via environment variables.
MAX_SEQS = int(os.environ.get("MFE_MAX_SEQS", "50"))
MAX_LEN = int(os.environ.get("MFE_MAX_LEN", "500"))

VIENNA = RNA.__version__

LAYOUTS = {
    "puzzler": ("Puzzler (default)", RNA.PLOT_TYPE_PUZZLER),
    "circular": ("Circular", RNA.PLOT_TYPE_CIRCULAR),
    "naview": ("NAView", RNA.PLOT_TYPE_NAVIEW),
}


def fold(seq, material, temp_c):
    """Return (structure, MFE kcal/mol, per-base paired probability list).

    The probabilities come from the partition function: for base i it is the
    sum over j of P(i pairs with j) in the thermodynamic ensemble (not only
    the MFE structure).
    """
    n = len(seq)
    with _lock:
        if material == "DNA":
            RNA.params_load_DNA_Mathews2004()
        else:
            RNA.params_load_RNA_Turner2004()
        md = RNA.md()
        md.temperature = temp_c
        # ViennaRNA expects U internally; T is just DNA notation for U.
        fc = RNA.fold_compound(seq.replace("T", "U"), md)
        structure, mfe = fc.mfe()
        fc.exp_params_rescale(mfe)
        fc.pf()
        bpp = fc.bpp()
    paired = [0.0] * n
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            p = bpp[i][j]
            paired[i - 1] += p
            paired[j - 1] += p
    paired = [min(1.0, max(0.0, p)) for p in paired]
    return structure, round(mfe, 2), paired


def pair_list(structure):
    stack, pairs = [], []
    for i, ch in enumerate(structure):
        if ch == "(":
            stack.append(i)
        elif ch == ")":
            pairs.append((stack.pop(), i))
    return pairs


def layout_xy(structure, layout):
    with _lock:
        RNA.cvar.rna_plot_type = LAYOUTS.get(layout, LAYOUTS["puzzler"])[1]
        c = RNA.get_xy_coordinates(structure)
        return [(c.get(i).X, c.get(i).Y) for i in range(len(structure))]


def shade(p):
    """Light yellow (p=0, mostly unpaired) -> dark blue (p=1, always paired)."""
    a, b = (255, 245, 214), (37, 99, 176)
    rgb = tuple(round(a[i] + (b[i] - a[i]) * p) for i in range(3))
    return "#%02x%02x%02x" % rgb, ("#ffffff" if p > 0.55 else "#111111")


def draw_structure(seq, structure, paired, layout):
    """Inline SVG with bases coloured by pairing probability.

    Built here from numbers and the letters A/C/G/T/U only (input is validated
    upstream), so it is safe to mark as Markup.
    """
    n = len(seq)
    xy = layout_xy(structure, layout)
    steps = sorted(math.dist(xy[i], xy[i + 1]) for i in range(n - 1)) or [1.0]
    k = 22.0 / (steps[len(steps) // 2] or 1.0)
    pts = [(x * k, -y * k) for x, y in xy]
    pad = 20
    minx = min(p[0] for p in pts) - pad
    miny = min(p[1] for p in pts) - pad
    w = max(p[0] for p in pts) + pad - minx
    h = max(p[1] for p in pts) + pad - miny

    parts = []
    for i in range(n - 1):
        (x1, y1), (x2, y2) = pts[i], pts[i + 1]
        parts.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#9aa" stroke-width="2"/>')
    for i, j in pair_list(structure):
        (x1, y1), (x2, y2) = pts[i], pts[j]
        parts.append(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="#d33" stroke-width="2.5"/>')
    for i, (x, y) in enumerate(pts):
        fill, ink = shade(paired[i])
        parts.append(
            f'<circle cx="{x:.1f}" cy="{y:.1f}" r="9" fill="{fill}" stroke="#555" stroke-width="0.8">'
            f'<title>{seq[i]}{i + 1}: P(paired) = {paired[i]:.2f}</title></circle>'
            f'<text x="{x:.1f}" y="{y + 3.8:.1f}" text-anchor="middle" font-size="11" '
            f'font-family="sans-serif" fill="{ink}" pointer-events="none">{seq[i]}</text>'
        )
    vb = f"{minx:.1f} {miny:.1f} {w:.1f} {h:.1f}"
    return Markup(f'<svg class="fold" viewBox="{vb}" xmlns="http://www.w3.org/2000/svg">' + "".join(parts) + "</svg>")


def parse_input(text):
    """Parse lines as 'SEQ' or 'name SEQ' / 'name,SEQ'. Returns (rows, errors)."""
    rows, errors = [], []
    for i, line in enumerate(text.splitlines(), 1):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = re.split(r"[\s,;]+", line)
        seq = parts[-1].upper()
        name = " ".join(parts[:-1]) or f"Seq {len(rows) + 1}"
        if not re.fullmatch(r"[ACGTU]+", seq):
            errors.append(f"Line {i}: only A, C, G, T, U allowed.")
            continue
        if len(seq) > MAX_LEN:
            errors.append(f"Line {i}: longer than {MAX_LEN} nt.")
            continue
        rows.append((name, seq))
    if len(rows) > MAX_SEQS:
        errors.append(f"At most {MAX_SEQS} sequences per request.")
        rows = rows[:MAX_SEQS]
    return rows, errors


def analyse(rows, material, temp_c, layout="puzzler"):
    out = []
    for name, seq in rows:
        if material == "DNA" and "U" in seq:
            seq = seq.replace("U", "T")
        elif material == "RNA" and "T" in seq:
            seq = seq.replace("T", "U")
        structure, mfe, paired = fold(seq, material, temp_c)
        pairs = structure.count("(")
        gc = (seq.count("G") + seq.count("C")) / len(seq)
        out.append(
            dict(name=name, seq=seq, length=len(seq), gc=gc, structure=structure,
                 mfe=mfe, pairs=pairs, paired_bases=2 * pairs,
                 svg=draw_structure(seq, structure, paired, layout))
        )
    return out


PAGE = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>MFE folding (ViennaRNA)</title>
<style>
 :root{--mono:ui-monospace,SFMono-Regular,Menlo,Consolas,monospace}
 body{font-family:system-ui,-apple-system,"Segoe UI",Roboto,sans-serif;font-size:15px;line-height:1.5;max-width:980px;margin:2rem auto;padding:0 1rem;color:#222}
 select,input,button,textarea{font:inherit;color:inherit}
 select,input[type=text],input:not([type]){padding:.25rem .4rem;border:1px solid #aaa;border-radius:4px;background:#fff}
 code{font-family:var(--mono);font-size:.9em;background:#f3f3f3;padding:0 .25em;border-radius:3px}
 .note{background:#f3f6fa;border:1px solid #d5dde8;padding:.6rem 1rem;border-radius:6px;color:#333}
 .err{background:#fde2e2;border:1px solid #d66;padding:.5rem 1rem;border-radius:6px;margin:.5rem 0}
 textarea{width:100%;height:9rem;font-family:var(--mono);font-size:.9rem;padding:.5rem;border:1px solid #aaa;border-radius:4px;box-sizing:border-box}
 table{border-collapse:collapse;width:100%;margin-top:1rem}
 th,td{border:1px solid #ccc;padding:.3rem .5rem;text-align:left}
 td.m{font-family:var(--mono);font-size:.85rem;word-break:break-all}
 tr.avg td{font-weight:600;background:#f3f3f3}
 label{margin-right:1rem} button{padding:.4rem 1.1rem;border:1px solid #888;border-radius:4px;background:#f4f4f4;cursor:pointer}
 button:hover{background:#e8e8e8}
 small{color:#555}
 .figs{display:flex;flex-wrap:wrap;gap:1rem;margin-top:1rem}
 figure{margin:0;border:1px solid #ddd;border-radius:6px;padding:.4rem;text-align:center;width:min(100%,440px)}
 svg.fold{width:100%;height:auto;max-height:440px;cursor:grab;touch-action:none;background:#fff}
 svg.fold:active{cursor:grabbing}
 figcaption{font-size:.85rem;color:#444}
 .legend{display:flex;align-items:center;gap:.5rem;font-size:.85rem;margin-top:1rem}
 .bar{width:160px;height:12px;border:1px solid #888;background:linear-gradient(90deg,#fff5d6,#2563b0)}
</style></head><body>
<h1>MFE folding tool</h1>
<p class="note">Single-strand MFE folding with ViennaRNA {{ vienna }}
(Turner 2004 parameters for RNA, Mathews 2004 for DNA). Free energies can differ from NUPACK's.</p>
<form method="post">
 <p><label>Material
  <select name="material">
   <option {{ 'selected' if material=='DNA' }}>DNA</option>
   <option {{ 'selected' if material=='RNA' }}>RNA</option>
   <option value="BOTH" {{ 'selected' if material=='BOTH' }}>Both (DNA and RNA)</option>
  </select></label>
  <label>Temperature (&deg;C)
  <input name="temp" value="{{ temp }}" size="5"></label>
  <label>Drawing layout
  <select name="layout">
   {% for key, (label, _) in layouts.items() %}
   <option value="{{ key }}" {{ 'selected' if layout==key }}>{{ label }}</option>
   {% endfor %}
  </select></label></p>
 <p>Sequences, one per line (<code>ACGU...</code> or <code>name ACGU...</code>). T and U are accepted either way.</p>
 <textarea name="seqs" spellcheck="false">{{ seqs }}</textarea>
 <p><button>Fold</button></p>
</form>
{% for e in errors %}<div class="err">{{ e }}</div>{% endfor %}
{% if results %}
<div class="legend"><span>Base colour = probability of being paired in the ensemble:</span>
 <span>0</span><span class="bar"></span><span>1</span>
 <small>&nbsp;Scroll to zoom, drag to pan, double-click to reset. Hover a base for its value.</small></div>
{% endif %}
{% for mat, rows in results %}
 <h2>{{ mat }}</h2>
 <table>
  <tr><th>Name</th><th>Length</th><th>GC</th><th>MFE structure (dot-bracket)</th>
      <th>&Delta;G (kcal/mol)</th><th>Base pairs</th><th>Paired bases</th></tr>
  {% for r in rows %}
  <tr><td>{{ r.name }}</td><td>{{ r.length }}</td><td>{{ '%.2f'|format(r.gc) }}</td>
      <td class="m">{{ r.seq }}<br>{{ r.structure }}</td>
      <td>{{ '%.2f'|format(r.mfe) }}</td><td>{{ r.pairs }}</td><td>{{ r.paired_bases }}</td></tr>
  {% endfor %}
  {% if rows|length > 1 %}
  <tr class="avg"><td colspan="4">Average over {{ rows|length }} sequences</td>
   <td>{{ '%.2f'|format(rows|sum(attribute='mfe')/rows|length) }}</td>
   <td>{{ '%.2f'|format(rows|sum(attribute='pairs')/rows|length) }}</td>
   <td>{{ '%.2f'|format(rows|sum(attribute='paired_bases')/rows|length) }}</td></tr>
  {% endif %}
 </table>
 <div class="figs">
  {% for r in rows %}
  <figure>{{ r.svg }}
   <figcaption>{{ r.name }} &middot; {{ mat }} &middot; {{ '%.2f'|format(r.mfe) }} kcal/mol</figcaption></figure>
  {% endfor %}
 </div>
{% endfor %}
{% if results %}<p><small>"Base pairs" = number of pairs (count of "(").
"Paired bases" = 2 &times; pairs. Check which one your exercise asks for.
Engine: ViennaRNA, Institute for Theoretical Chemistry, University of Vienna.</small></p>{% endif %}
<script>
document.querySelectorAll('svg.fold').forEach(function (svg) {
  var o = svg.getAttribute('viewBox').split(' ').map(Number), v = o.slice(), start = null;
  function set() { svg.setAttribute('viewBox', v.join(' ')); }
  svg.addEventListener('wheel', function (e) {
    e.preventDefault();
    var r = svg.getBoundingClientRect(), px = (e.clientX - r.left) / r.width, py = (e.clientY - r.top) / r.height;
    var f = e.deltaY < 0 ? 0.85 : 1 / 0.85, nw = v[2] * f, nh = v[3] * f;
    if (nw >= o[2]) { v = o.slice(); set(); return; }
    if (nw < o[2] / 25) return;
    v = [v[0] + (v[2] - nw) * px, v[1] + (v[3] - nh) * py, nw, nh]; set();
  }, { passive: false });
  svg.addEventListener('pointerdown', function (e) { start = { x: e.clientX, y: e.clientY, v: v.slice() }; svg.setPointerCapture(e.pointerId); });
  svg.addEventListener('pointermove', function (e) {
    if (!start) return;
    var r = svg.getBoundingClientRect();
    v = [start.v[0] - (e.clientX - start.x) * start.v[2] / r.width, start.v[1] - (e.clientY - start.y) * start.v[3] / r.height, start.v[2], start.v[3]];
    set();
  });
  svg.addEventListener('pointerup', function () { start = null; });
  svg.addEventListener('dblclick', function () { v = o.slice(); set(); });
});
</script>
</body></html>"""


@app.route("/", methods=["GET", "POST"])
def index():
    material = request.form.get("material", "RNA")
    seqs = request.form.get("seqs", "")
    temp_raw = request.form.get("temp", "37")
    layout = request.form.get("layout", "puzzler")
    if layout not in LAYOUTS:
        layout = "puzzler"
    results, errors = [], []
    if request.method == "POST":
        try:
            temp = float(temp_raw)
            if not -10 <= temp <= 100:
                raise ValueError
        except ValueError:
            temp, temp_raw = 37.0, "37"
            errors.append("Temperature must be a number between -10 and 100; using 37.")
        rows, perr = parse_input(seqs)
        errors += perr
        if not rows and not perr:
            errors.append("Enter at least one sequence.")
        mats = ["DNA", "RNA"] if material == "BOTH" else [material if material in ("DNA", "RNA") else "RNA"]
        if rows:
            results = [(m, analyse(rows, m, temp, layout)) for m in mats]
    return render_template_string(PAGE, material=material, seqs=seqs, temp=temp_raw,
                                  layout=layout, layouts=LAYOUTS, results=results, errors=errors,
                                  vienna=VIENNA)


if __name__ == "__main__":
    # Local development only (localhost). Hosted deployments use gunicorn, see Dockerfile.
    app.run(host="127.0.0.1", port=5050)
