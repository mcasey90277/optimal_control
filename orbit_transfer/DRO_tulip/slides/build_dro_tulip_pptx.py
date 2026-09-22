"""Build the 5-slide DRO -> tulip overview deck (black, 16:9).

Run:  ~/ai_council/venv/bin/python build_dro_tulip_pptx.py
Needs assets/ from make_slide_assets.m (movie, CR3BP figure) and the
library-of-record phase-torus figure. Rebuilding overwrites manual edits.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.text import PP_ALIGN

HERE = os.path.dirname(os.path.abspath(__file__))
A = os.path.join(HERE, "assets")
LIB = os.path.join(HERE, "..", "indirect", "results", "library_70mN_24x24_final")

BLACK = RGBColor(0, 0, 0)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
GOLD = RGBColor(0xFF, 0xD8, 0x4D)
BLUE = RGBColor(0x6C, 0xB8, 0xFF)
GREY = RGBColor(0xB0, 0xB0, 0xB0)
GREEN = RGBColor(0x5C, 0xD6, 0x7A)

# numbers from the render (assets/extremes_70mN.txt)
EXT = {}
with open(os.path.join(A, "extremes_70mN.txt")) as f:
    for line in f:
        tok = line.split()
        d = dict(t.split("=") for t in tok[1:] if "=" in t)
        EXT[tok[0]] = d


def eq_png(name, lines, fs=22, w=6.0, h=2.0):
    """Render mathtext lines to a white-on-black PNG."""
    path = os.path.join(A, name)
    fig = plt.figure(figsize=(w, h), facecolor="black")
    n = len(lines)
    for k, s in enumerate(lines):
        fig.text(0.02, 1 - (k + 0.6) / n, s, color="white", fontsize=fs,
                 va="center", ha="left")
    fig.savefig(path, dpi=200, facecolor="black")
    plt.close(fig)
    return path


def new_slide(prs, title, sub=None):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    s.background.fill.solid()
    s.background.fill.fore_color.rgb = BLACK
    tb = s.shapes.add_textbox(Inches(0.5), Inches(0.25), Inches(12.3), Inches(0.8))
    tf = tb.text_frame
    tf.word_wrap = True
    r = tf.paragraphs[0].add_run()
    r.text = title
    r.font.size = Pt(30)
    r.font.bold = True
    r.font.color.rgb = WHITE
    if sub:
        p = tf.add_paragraph()
        r = p.add_run()
        r.text = sub
        r.font.size = Pt(17)
        r.font.color.rgb = GREY
    return s


def bullets(s, x, y, w, h, items, size=18):
    """items: list of (level, [(text, color, bold), ...]) or plain strings."""
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    first = True
    for it in items:
        if isinstance(it, str):
            it = (0, [(it, WHITE, False)])
        lvl, runs = it
        p = tf.paragraphs[0] if first else tf.add_paragraph()
        first = False
        p.space_after = Pt(8)
        lead = p.add_run()
        lead.text = ("• " if lvl == 0 else "     – ")
        lead.font.size = Pt(size - 2 * lvl)
        lead.font.color.rgb = GOLD if lvl == 0 else GREY
        for text, col, bold in runs:
            r = p.add_run()
            r.text = text
            r.font.size = Pt(size - 2 * lvl)
            r.font.color.rgb = col
            r.font.bold = bold
    return tb


def table(s, x, y, w, rows, colw, size=15, rowh=0.42):
    t = s.shapes.add_table(len(rows), len(rows[0]), Inches(x), Inches(y),
                           Inches(w), Inches(rowh * len(rows))).table
    for c, cw in enumerate(colw):
        t.columns[c].width = Inches(cw)
    for r_, row in enumerate(rows):
        for c, val in enumerate(row):
            cell = t.cell(r_, c)
            cell.fill.solid()
            cell.fill.fore_color.rgb = RGBColor(0x22, 0x22, 0x22) if r_ == 0 else BLACK
            p = cell.text_frame.paragraphs[0]
            p.text = ""
            run = p.add_run()
            run.text = val
            run.font.size = Pt(size)
            run.font.bold = (r_ == 0 or c == 0)
            run.font.color.rgb = GOLD if r_ == 0 else (WHITE if c == 0 else BLUE)
    return t


prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)

fast, slow = EXT["FASTEST"], EXT["SLOWEST"]

# ---------------------------------------------------------------- slide 1
s = new_slide(prs, "Minimum-time transfer: thrust all the way, point along the primer",
              "A 70 mN electric thruster moves a 150 kg spacecraft from a lunar DRO to a 7-petal tulip in 16-22 days")
eq1 = eq_png("eq_mintime.png", [
    r"$\min\ t_f \quad \mathrm{s.t.}\quad \dot r = v,\ \ \dot v = g_{\mathrm{CR3BP}}(r,v) + \frac{T}{m}\,\hat u,\ \ \dot m = -\frac{T}{c}$",
    r"$\hat u^{\ast} = -\lambda_v/\|\lambda_v\|,\qquad \|u\| = 1\ \ \mathrm{(all\ burn)}$",
    r"$H(t_f) = 0,\qquad \mathrm{unknowns}\ \ z_8 = [\lambda_r;\ \lambda_v;\ \lambda_m;\ t_f]$",
], fs=19, w=9.6, h=2.1)
s.shapes.add_picture(eq1, Inches(0.5), Inches(1.45), Inches(7.2))
table(s, 8.1, 1.5, 4.8, [
    ["Engine / spacecraft", "Value"],
    ["Thrust T", "70 mN"],
    ["Specific impulse", "900 s  (c = 8.83 km/s)"],
    ["Initial mass m₀", "150 kg"],
    ["Acceleration T/m₀", "0.47 mm/s²"],
    ["Propellant flow", "0.69 kg/day"],
], colw=[2.3, 2.5], size=15)
bullets(s, 0.5, 3.55, 7.3, 3.8, [
    (0, [("Pontryagin: ", GOLD, True), ("the control is fixed by the costates; the problem becomes an 8-unknown two-point BVP", WHITE, False)]),
    (0, [("Solved in two stages: ", GOLD, True), ("direct collocation (CasADi/IPOPT) finds the basin; multiple shooting on the PMP field polishes the costates to ~1e-13", WHITE, False)]),
    (0, [("Independent witness: ", GOLD, True), ("pumpkyn's own tfMin accepts every stored z₈ unchanged", WHITE, False)]),
    (0, [("Fuel follows time: ", GOLD, True), ("all-burn (a theorem for min time) means propellant = 0.69 kg/day × t_f, so across the library the fastest transfer is also the cheapest", WHITE, False)]),
], size=17)
s.shapes.add_picture(os.path.join(HERE, "..", "indirect", "results", "transfer_3d_anchor.png"),
                     Inches(8.1), Inches(4.3), Inches(4.8))

# ---------------------------------------------------------------- slide 2
s = new_slide(prs, "CR3BP: both orbits are periodic in the Earth-Moon rotating frame",
              "Circular restricted three-body problem, nondimensionalised on the Earth-Moon distance and period")
eq2 = eq_png("eq_cr3bp.png", [
    r"$\ddot x - 2\dot y = \partial\Omega/\partial x,\quad \ddot y + 2\dot x = \partial\Omega/\partial y,\quad \ddot z = \partial\Omega/\partial z$",
    r"$\Omega = \frac{1}{2}(x^2+y^2) + \frac{1-\mu}{r_1} + \frac{\mu}{r_2},\qquad C = 2\Omega - \|v\|^2$",
], fs=19, w=7.6, h=1.3)
s.shapes.add_picture(eq2, Inches(0.5), Inches(1.45), Inches(6.4))
table(s, 0.5, 2.75, 6.3, [
    ["Constant", "Value"],
    ["Mass ratio μ", "0.0121506"],
    ["Length unit l*", "389,703 km"],
    ["Time unit t*", "382,981 s = 4.43 d"],
    ["Departure: DRO", "period τ = 1 (4.43 d)"],
    ["Arrival: 7-petal tulip", "period 5.236 (23.2 d)"],
], colw=[2.9, 3.4], size=15)
bullets(s, 0.5, 5.45, 6.4, 2.0, [
    (0, [("DRO: ", GOLD, True), ("stable, retrograde, planar, around the Moon", WHITE, False)]),
    (0, [("Tulip: ", GOLD, True), ("7-petal 3D resonant orbit over the lunar poles", WHITE, False)]),
    (0, [("Orbits from pumpkyn's catalogued families; ", GOLD, True), ("positions set by phase fractions s_D, s_A ∈ [0,1)", WHITE, False)]),
], size=16)
s.shapes.add_picture(os.path.join(A, "cr3bp_geometry.png"), Inches(7.0), Inches(1.45), Inches(6.0))

# ---------------------------------------------------------------- slide 3
s = new_slide(prs, "The costate library: every phase pair solved, certified and audited",
              "A lookup table of converged PMP costates, so any transfer starts from a known root instead of a cold guess")
bullets(s, 0.5, 1.5, 7.0, 5.8, [
    (0, [("What an entry is: ", GOLD, True), ("(s_D, s_A) → z₈ = [λ(7); t_f], plus its certificate", WHITE, False)]),
    (0, [("70 mN library of record: ", GOLD, True), ("24 × 24 phase grid, ", WHITE, False), ("576 of 576 cells certified", GREEN, True)]),
    (1, [("rebuilt unattended by one script in 24 h; fail-closed audit 576 OK / 0 bad", WHITE, False)]),
    (1, [("five solution families found; the library keeps the fastest root in each cell", WHITE, False)]),
    (0, [("Necessary conditions: ", GOLD, True), ("flown miss, H = 0, transversality, adjoint equations, minimum-principle gap", WHITE, False)]),
    (0, [("Sufficiency: ", GOLD, True), ("conjugate-point test (no interior crossing), Legendre / switching / H6 margins, abnormal-lift rank", WHITE, False)]),
    (0, [("Wider programme: ", GOLD, True), ("~18,400 min-time entries over DRO / halo / DPO → tulip and halo ↔ halo, thrust 0.5-15 N", WHITE, False)]),
    (0, [("Now: ", GOLD, True), ("interpolating between cells. Along arrival phase at 1/96 spacing, 87% of blended guesses converge (70% from the nearest entry); departure spacing still to be measured", WHITE, False)]),
], size=17)
table(s, 7.9, 1.6, 5.0, [
    ["Pipeline stage", "Tool"],
    ["1  Basin", "direct collocation (IPOPT)"],
    ["2  Costates", "dual harvest → multiple shooting"],
    ["3  Continuation", "pseudo-arclength along phase"],
    ["4  Certify", "certify_root gate stack"],
    ["5  Witness", "pumpkyn tfMin"],
    ["6  Audit", "re-fly every entry from its keys"],
], colw=[2.1, 2.9], size=14)

# ---------------------------------------------------------------- slide 4
s = new_slide(prs, f"Phasing alone costs {float(slow['tf']) - float(fast['tf']):.1f} days: "
                   f"fastest {float(fast['tf']):.2f} d vs slowest {float(slow['tf']):.2f} d",
              "Same engine, same two orbits, only the departure/arrival phases differ; both panels run on one clock")
gif = os.path.join(A, "extremes_70mN.gif")
GW = 10.4
s.shapes.add_picture(gif, Inches((13.333 - GW) / 2), Inches(1.4), Inches(GW), Inches(GW * 720 / 1280))

# ---------------------------------------------------------------- slide 5
s = new_slide(prs, "Arrival phase sets the transfer time; departure phase mostly does not",
              "Minimum time over the 24 × 24 phase torus (70 mN, Isp 900 s, 150 kg)")
s.shapes.add_picture(os.path.join(LIB, "phase_torus_70mN.png"), Inches(0.4), Inches(1.35), Inches(6.9))
bullets(s, 7.6, 1.6, 5.4, 5.6, [
    (0, [("Vertical stripes: ", GOLD, True), ("t_f is almost constant down a column, so where you meet the tulip dominates", WHITE, False)]),
    (0, [("Range: ", GOLD, True), (f"{float(fast['tf']):.2f} to {float(slow['tf']):.2f} d over the torus (about 33%)", WHITE, False)]),
    (0, [("Bottom row and blocks: ", GOLD, True), ("sharp steps are hand-overs between solution families, not noise", WHITE, False)]),
    (0, [("Every cell: ", GOLD, True), ("certified root, second-order sweep clean (0 interior conjugate crossings)", WHITE, False)]),
    (0, [("Use: ", GOLD, True), ("pick the arrival phase for schedule; the library returns the costates to fly it", WHITE, False)]),
], size=17)

out = os.path.join(HERE, "DRO_tulip_overview.pptx")
prs.save(out)
for k, sl in enumerate(prs.slides, 1):
    kinds = sorted({str(sh.shape_type) for sh in sl.shapes})
    print(k, kinds)
print("saved", out)
