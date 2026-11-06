import re
from app import app, analyse, parse_input

SEQS = """AGTTGTGACTTTTAGGAGCCCCCCCTAAAGAGGAGCTTAC
ACCCGGGTTTACTGCCTCTTCGAAGTACGCGAAAGTTAAG
AGTACTTTCGCTGAAATCCGGTAGGACATTCGGGACACTC
TGCAAGTATTCCGCAGGGTTATAGCACACTGGCGTTACAC
TCACCAGACGATCCATTTCAGAGCTGTCGTTGAAGTCGAG"""

rows, errs = parse_input(SEQS)
assert len(rows) == 5 and not errs
for mat in ("DNA", "RNA"):
    res = analyse(rows, mat, 37.0)
    for r in res:
        assert len(r["structure"]) == 40 and r["gc"] == 0.5
        assert r["structure"].count("(") == r["structure"].count(")")
        assert r["mfe"] <= 0
    n = len(res)
    print(mat, "mfe:", [r["mfe"] for r in res], "avg", round(sum(r["mfe"] for r in res)/n, 2),
          "pairs:", [r["pairs"] for r in res], "avg", sum(r["pairs"] for r in res)/n)

# T/U handling: same sequence given as T or U must give identical RNA result
a = analyse([("x", "GGGAAACCC")], "RNA", 37)[0]; b = analyse([("x", "GGGAAACCC".replace("T","U"))], "RNA", 37)[0]
assert a["mfe"] == b["mfe"]
# temperature matters
assert analyse(rows[:1], "RNA", 60)[0]["mfe"] > analyse(rows[:1], "RNA", 37)[0]["mfe"]
# bad input
_, e = parse_input("ACGX\nAC GT")
assert len(e) == 1, e

c = app.test_client()
assert c.get("/").status_code == 200
r = c.post("/", data=dict(material="BOTH", temp="37", seqs=SEQS))
assert r.status_code == 200 and b"Average over 5" in r.data and b"<h2>DNA</h2>" in r.data
r = c.post("/", data=dict(material="RNA", temp="abc", seqs="GGGAAACCC"))
assert b"Temperature must be" in r.data
r = c.post("/", data=dict(material="RNA", temp="37", seqs=""))
assert b"at least one sequence" in r.data
print("ALL OK")

# structure drawings
from app import fold, LAYOUTS
for lay in LAYOUTS:
    res = analyse(rows, "RNA", 37.0, lay)
    for r in res:
        svg = str(r["svg"])
        assert svg.startswith("<svg") and svg.count("<circle") == 40, (lay, svg.count("<circle"))
        assert svg.count('stroke="#d33"') == r["pairs"]
# probabilities are valid and high for MFE-paired bases on a strong hairpin
st, e, pp = fold("GGGGGAAAACCCCC", "RNA", 37.0)
assert len(pp) == 14 and all(0 <= x <= 1 for x in pp) and pp[0] > 0.9 and pp[7] < 0.1, pp
r = c.post("/", data=dict(material="BOTH", temp="37", seqs=SEQS, layout="circular"))
assert r.data.count(b"<figure>") == 10 and b'class="fold"' in r.data
r = c.post("/", data=dict(material="RNA", temp="37", seqs="GGGAAACCC", layout="evil"))
assert r.status_code == 200 and r.data.count(b"<figure>") == 1
print("DRAWINGS OK")

# input handling
_, e = parse_input("ACGX")
assert len(e) == 1 and "X" in e[0], e
_, e = parse_input(">seq1")
assert len(e) == 1 and "FASTA" in e[0], e
_, e = parse_input("\n".join(["x"] * 30))
assert len(e) == 11 and "20 more" in e[-1], e          # 10 listed + summary line
r = c.post("/", data=dict(material="RNA", temp="37", seqs="x" * (300 * 1024)))
assert r.status_code == 413 and b"too large" in r.data   # oversized body rejected up front
r = c.post("/", data=dict(material="RNA", temp="37", seqs="<b>x</b> GGGAAACCC"))
assert b"<b>x</b>" not in r.data and b"&lt;b&gt;x&lt;/b&gt;" in r.data   # names are escaped
assert b"placeholder=" in c.get("/").data
print("INPUT OK")
