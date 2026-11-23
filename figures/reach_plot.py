#!/usr/bin/env python3
"""Figure: where the tests reach, over all 2320 Smith parameter sets with r <= 135.
Reads ../smith_complete_r135.tsv (from scan_complete.py; its 'checks' column is the single
feasibility predicate used in the paper) and writes reach.pdf / reach.png.
Usage (from this directory): MPLCONFIGDIR=<writable dir> python -B reach_plot.py

y = E/(r(r+1)) - 1 = 2u/(1-2u), so a fixed offset u is a horizontal line.
Curves: Kre(t) has y = 2/(r+1); series (i) has y = 4/r; GQ(q,q^2) has y = 2/r.
"""
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE = "#2a78d6", "#eb6834"          # reference palette slots 1, 2
INK, MUTED, FAINT = "#2b2b29", "#8a897f", "#d9d8d0"

rows = list(csv.DictReader(open("../smith_complete_r135.tsv"), delimiter="\t"))


def category(R):
    ex = "excl" in (R["A"], R["V"])
    if ex and R["checks"] == "True":
        return "second" if "True" in (R["S2P"], R["S2S"]) else "new"
    if ex:
        return "infeasible"
    if "pass" in (R["A"], R["V"]):
        return "pass"
    if "undet" in (R["A"], R["V"]):
        return "undet"
    return "na"


style = {  # drawn in this order (background first)
    "na":         dict(s=4, c=FAINT, marker="o", lw=0,
                       label="outside both applicability conditions"),
    "undet":      dict(s=6, facecolors="none", edgecolors=FAINT, marker="o", lw=0.5,
                       label=r"$\Delta$ not determined by the parameters"),
    "pass":       dict(s=9, c=MUTED, marker="o", lw=0, label="congruence satisfied"),
    "infeasible": dict(s=18, c=MUTED, marker="x", lw=0.8,
                       label="excluded; fails the feasibility checks"),
    "second":     dict(s=26, c=ORANGE, marker="^", lw=0.4, edgecolors="white",
                       label="excluded; also by Prop. 4.4"),
    "new":        dict(s=30, c=BLUE, marker="o", lw=0.5, edgecolors="white",
                       label="excluded; no Prop. 4.4 exclusion"),
}


def ycoord(r, s):
    return (r - s) / (r * (r + 1)) - 1          # = 2u/(1-2u)


pts = {k: ([], []) for k in style}
real = ([], [])
headline = {}
for R in rows:
    n, k, r, s = int(R["n"]), int(R["k"]), int(R["r"]), int(R["s"])
    y = ycoord(r, s)
    pts[category(R)][0].append(r)
    pts[category(R)][1].append(y)
    if R["realized"] == "True":
        real[0].append(r)
        real[1].append(y)
    if n in (392, 784, 1080):
        headline[n] = (r, y)

fig, ax = plt.subplots(figsize=(6.4, 4.3), dpi=150)
for key, st in style.items():
    xs, ys = pts[key]
    ax.scatter(xs, ys, zorder=3 if key in ("new", "second") else 2, **st)
ax.scatter(*real, s=70, marker="*", facecolors="none", edgecolors=INK, lw=0.7, zorder=4,
           label="realized in the table")
rs = [x / 10 for x in range(10, 1360)]
ax.plot(rs, [2 / (x + 1) for x in rs], color=MUTED, lw=0.6, ls="--", zorder=1)
ax.text(105, 2 / 106 * 0.78, r"$\mathrm{Kre}(t)$", fontsize=7, color=MUTED, ha="center", va="top")
labels = {392: "392", 784: r"784 = Kre(4)", 1080: "1080"}
offsets = {392: (-30, 10), 784: (14, -18), 1080: (10, 8)}
for n, (r, y) in headline.items():
    ax.annotate(labels[n], xy=(r, y), xytext=offsets[n], textcoords="offset points",
                fontsize=7, color=INK, arrowprops=dict(arrowstyle="-", color=INK, lw=0.5))
ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlabel(r"positive eigenvalue $r$", fontsize=8, color=INK)
ax.set_ylabel(r"$2u/(1-2u)$", fontsize=8, color=INK)
for sp in ("top", "right"):
    ax.spines[sp].set_visible(False)
for sp in ("left", "bottom"):
    ax.spines[sp].set_color(MUTED)
ax.tick_params(colors=MUTED, labelcolor=INK, labelsize=7)
leg = fig.legend(fontsize=6.5, frameon=False, loc="lower center", ncol=3, handletextpad=0.3,
                 columnspacing=1.2)
for t in leg.get_texts():
    t.set_color(INK)
fig.tight_layout(rect=(0, 0.1, 1, 1))
fig.savefig("reach.pdf")
fig.savefig("reach.png", dpi=150)
print({k: len(v[0]) for k, v in pts.items()}, "realized:", len(real[0]), "headline:", headline)
