"""Build the DRO -> tulip overview deck (black, 16:9): title + 5 content
slides (v3: Mike's v2 order and wording, 2026-09-22), and an optional
4-slide appendix on the optimality checks (APPENDIX flag; v2 dropped it).

Run:  ~/ai_council/venv/bin/python build_dro_tulip_pptx.py
Needs assets/ from make_slide_assets.m (movie, CR3BP figure, poster frame)
and the library-of-record phase-torus figure. Rebuilding overwrites manual
edits: Mike's 2026-09-22 edits are baked in below; his edited files are
kept as assets/DRO_tulip_overview_user_2026-09-22.pptx and
assets/DRO_tulip_overview_v2_user_2026-09-22.pptx.
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR

APPENDIX = False            # v2 removed the appendix; True restores it

HERE = os.path.dirname(os.path.abspath(__file__))
A = os.path.join(HERE, "assets")
LIB = os.path.join(HERE, "..", "indirect", "results", "library_70mN_24x24_final")

BLACK = RGBColor(0, 0, 0)
WHITE = RGBColor(0xFF, 0xFF, 0xFF)
GOLD = RGBColor(0xFF, 0xD8, 0x4D)
BLUE = RGBColor(0x6C, 0xB8, 0xFF)
GREY = RGBColor(0xB0, 0xB0, 0xB0)
GREEN = RGBColor(0x5C, 0xD6, 0x7A)
DARK = RGBColor(0x22, 0x22, 0x22)

# numbers from the render (assets/extremes_70mN.txt)
EXT = {}
with open(os.path.join(A, "extremes_70mN.txt")) as f:
    for line in f:
        tok = line.split()
        EXT[tok[0]] = dict(t.split("=") for t in tok[1:] if "=" in t)
fast, slow = EXT["FASTEST"], EXT["SLOWEST"]


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


def black_slide(prs):
    s = prs.slides.add_slide(prs.slide_layouts[6])
    s.background.fill.solid()
    s.background.fill.fore_color.rgb = BLACK
    return s


def new_slide(prs, title, sub=None, tag=None):
    s = black_slide(prs)
    tb = s.shapes.add_textbox(Inches(0.5), Inches(0.25), Inches(12.3), Inches(0.85))
    tf = tb.text_frame
    tf.word_wrap = True
    r = tf.paragraphs[0].add_run()
    r.text = title
    r.font.size = Pt(30)
    r.font.bold = True
    r.font.color.rgb = WHITE
    if sub:
        r = tf.add_paragraph().add_run()
        r.text = sub
        r.font.size = Pt(17)
        r.font.color.rgb = GREY
    if tag:                                   # appendix marker, top right
        t = s.shapes.add_textbox(Inches(11.6), Inches(7.0), Inches(1.6), Inches(0.4))
        r = t.text_frame.paragraphs[0].add_run()
        r.text = tag
        r.font.size = Pt(12)
        r.font.color.rgb = GREY
    return s


def bullets(s, x, y, w, h, items, size=18):
    """items: list of (level, [(text, color, bold), ...])."""
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    for k, (lvl, runs) in enumerate(items):
        p = tf.paragraphs[0] if k == 0 else tf.add_paragraph()
        p.space_after = Pt(8 if lvl == 0 else 4)
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


def table(s, x, y, w, rows, colw, size=15, rowh=0.42, colors=None):
    t = s.shapes.add_table(len(rows), len(rows[0]), Inches(x), Inches(y),
                           Inches(w), Inches(rowh * len(rows))).table
    for c, cw in enumerate(colw):
        t.columns[c].width = Inches(cw)
    colors = colors or [WHITE] + [BLUE] * (len(rows[0]) - 1)
    for r_, row in enumerate(rows):
        for c, val in enumerate(row):
            cell = t.cell(r_, c)
            cell.fill.solid()
            cell.fill.fore_color.rgb = DARK if r_ == 0 else BLACK
            cell.margin_top = cell.margin_bottom = Inches(0.04)
            p = cell.text_frame.paragraphs[0]
            p.text = ""
            run = p.add_run()
            run.text = val
            run.font.size = Pt(size)
            run.font.bold = (r_ == 0 or c == 0)
            run.font.color.rgb = GOLD if r_ == 0 else colors[c]
    return t


def rule(s, x, y, w, col=GOLD, h=0.05):
    b = s.shapes.add_shape(MSO_SHAPE.RECTANGLE, Inches(x), Inches(y), Inches(w), Inches(h))
    b.fill.solid()
    b.fill.fore_color.rgb = col
    b.line.fill.background()
    return b


def text(s, x, y, w, h, txt, size, col=WHITE, bold=False):
    tb = s.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tb.text_frame.word_wrap = True
    for k, line in enumerate(txt.split("\n")):
        p = tb.text_frame.paragraphs[0] if k == 0 else tb.text_frame.add_paragraph()
        r = p.add_run()
        r.text = line
        r.font.size = Pt(size)
        r.font.bold = bold
        r.font.color.rgb = col
    return tb


prs = Presentation()
prs.slide_width = Inches(13.333)
prs.slide_height = Inches(7.5)

# ================================================================ title
s = black_slide(prs)
s.shapes.add_picture(os.path.join(A, "title_scene.png"), Inches(6.35), Inches(1.05), Inches(6.8))
text(s, 0.7, 1.0, 5.6, 2.2, "DRO-to-Tulip\nMinimum-Time\nOrbit Transfers", 40, WHITE, True)
rule(s, 0.75, 3.55, 2.2)
text(s, 0.7, 3.75, 5.7, 0.8,
     "A costate library for low-thrust transfers in the Earth-Moon CR3BP", 20, GREY)
text(s, 0.7, 4.95, 5.7, 0.5, "70 mN  \u00b7  Isp 900 s  \u00b7  150 kg  \u00b7  576 phase pairs",
     16, GOLD, True)
text(s, 0.7, 6.65, 5.7, 0.4, "September 2026", 14, GREY)
text(s, 7.2, 6.75, 5.9, 0.4,
     f"Fastest certified transfer: {float(fast['tf']):.2f} d (blue), DRO (green) \u2192 tulip (red)",
     11, GREY)

# ================================================================ overview (movie)
s = new_slide(prs, "We Computed A Costate Library For DRO-to-Tulip Orbit Transfers")
GW = 9.0
s.shapes.add_picture(os.path.join(A, "extremes_70mN.gif"), Inches((13.333 - GW) / 2),
                     Inches(1.0), Inches(GW), Inches(GW * 720 / 1280))
dd = float(slow['tf']) - float(fast['tf'])
bullets(s, (13.333 - GW) / 2, 6.15, GW, 1.3, [
    (0, [("Costates computed on a grid of departure \u00d7 arrival phases: ", WHITE, False),
         ("24 \u00d7 24 = 576 phase pairs", GOLD, True)]),
    (0, [("Phasing alone costs ", WHITE, False), (f"{dd:.1f} days", GOLD, True),
         (f": fastest {float(fast['tf']):.2f} d (left) vs slowest {float(slow['tf']):.2f} d (right)", WHITE, False)]),
], size=17)

# ================================================================ library
s = new_slide(prs, "The costate library: costates computed for pairs of departure and arrival phases")
bullets(s, 0.5, 1.5, 7.1, 5.9, [
    (0, [("What an entry is: ", GOLD, True), ("(s_D, s_A) \u2192 z\u2088 = [\u03bb(7); t_f], plus its certificate", WHITE, False)]),
    (0, [("A root: ", GOLD, True), ("a z\u2088 that solves the shooting equations: flown from the DRO state with u = \u2212\u03bb_v/|\u03bb_v|, it arrives on the tulip state at t_f with H = 0 and \u03bb_m(t_f) = 0. Each root is one extremal", WHITE, False)]),
    (1, [("the equations are nonlinear, so a cell can hold several roots (distinct extremals, different t_f); goal is to find a local minimizer", WHITE, False)]),
    (0, [("70 mN costate library: ", GOLD, True), ("24 \u00d7 24 phase grid, ", WHITE, False), ("576 of 576 cells certified", GREEN, True)]),
    (1, [("five solution families found; the library keeps the fastest root in each cell", WHITE, False)]),
    (0, [("Necessary conditions: ", GOLD, True), ("flown miss, H = 0, transversality, adjoint equations, minimum-principle gap", WHITE, False)]),
    (0, [("Sufficiency: ", GOLD, True), ("conjugate-point test (no interior crossing), Legendre / switching / H6 margins, abnormal-lift rank", WHITE, False)]),
    (0, [("Wider effort underway: ", GOLD, True), ("~18,400 min-time entries over DRO / halo / DPO \u2192 tulip and halo \u2194 halo, thrust 0.5-15 N", WHITE, False)]),
    (0, [("Current investigation: ", GOLD, True), ("interpolating between cells. Along arrival phase at 1/96 spacing, 87% of blended guesses converge (70% from the nearest entry); departure spacing still to be measured", WHITE, False)]),
], size=16)
# the build as a five-step flow (counts from the catalog's entry notes,
# assets/provenance_counts.txt written by make_provenance_torus.m)
with open(os.path.join(A, "provenance_counts.txt")) as f:
    PC = {k: int(v) for k, v in (t.split("=") for t in f.read().split())}
FX, FW, FY, FH, GAP = 8.25, 4.6, 1.95, 0.82, 0.27
text(s, FX, FY - 0.5, FW, 0.4, "How the 576 entries were built", 17, GOLD, True)
steps = [("1  Seed roots", f"{PC['seed']} on the s_D = 0 row: Darin's root + direct solves"),
         ("2  Arclength along s_A", f"{PC['arc']} more cells on the s_D = 0 row"),
         ("3  Ribs along s_D", f"{PC['rib']} cells, multiple-shooting corrector each step"),
         ("4  Direct fills", f"{PC['direct']} cells where a rib stalled"),
         ("5  Certify + audit", "every one of the 576 cells")]
for k, (hd, sub) in enumerate(steps):
    y = FY + k*(FH + GAP)
    b = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, Inches(FX), Inches(y), Inches(FW), Inches(FH))
    b.adjustments[0] = 0.18
    b.fill.solid();  b.fill.fore_color.rgb = RGBColor(0x14, 0x1E, 0x2E)
    b.line.color.rgb = GOLD if k == len(steps) - 1 else BLUE;  b.line.width = Pt(1.75)
    tf = b.text_frame;  tf.word_wrap = True;  tf.vertical_anchor = MSO_ANCHOR.MIDDLE
    tf.margin_top = tf.margin_bottom = Inches(0.02)
    p = tf.paragraphs[0];  p.alignment = PP_ALIGN.CENTER
    r = p.add_run();  r.text = hd;  r.font.size = Pt(16);  r.font.bold = True;  r.font.color.rgb = WHITE
    p = tf.add_paragraph();  p.alignment = PP_ALIGN.CENTER
    r = p.add_run();  r.text = sub;  r.font.size = Pt(12);  r.font.color.rgb = GREY
    if k < len(steps) - 1:
        a = s.shapes.add_shape(MSO_SHAPE.DOWN_ARROW, Inches(FX + FW/2 - 0.15), Inches(y + FH + 0.03),
                               Inches(0.30), Inches(GAP - 0.06))
        a.fill.solid();  a.fill.fore_color.rgb = GOLD;  a.line.fill.background()

# ================================================================ provenance
nCont = PC['arc'] + PC['rib']
s = new_slide(prs, f"Six seed roots grow the library: {100*nCont/576:.0f}% of cells come from continuation",
              "Where every cell of the 24 \u00d7 24 library came from, read from each entry's own provenance note")
s.shapes.add_picture(os.path.join(A, "phase_torus_provenance.png"), Inches(0.35), Inches(1.35), Inches(7.7))
bullets(s, 8.35, 1.55, 4.75, 5.8, [
    (0, [(f"Seeds ({PC['seed']}): ", GOLD, True), ("all on the s_D = 0 row. One is Darin's pumpkynPie root, one came from an earlier phase sweep, four are direct solves", WHITE, False)]),
    (0, [(f"Arclength along s_A ({PC['arc']}): ", GOLD, True), ("pseudo-arclength walks from the seeds fill the rest of that row, through folds", WHITE, False)]),
    (0, [(f"Ribs along s_D ({PC['rib']}): ", GOLD, True), ("from each s_D = 0 cell, 23 steps in \u2212s_D (through the wrap at s_D = 1), so rows are shown in walk order", WHITE, False)]),
    (0, [(f"Direct fills ({PC['direct']}): ", GOLD, True), ("where a rib stalled, a direct solve warm-started from a certified neighbour; two columns are almost all fills", WHITE, False)]),
    (0, [("Same certificate for all: ", GOLD, True), ("provenance changes how a root was found, not how it is checked", WHITE, False)]),
], size=15)

# ================================================================ phase grid
s = new_slide(prs, "Transfer time more sensitive to arrival phase than departure phase",
              "Minimum time over the 24 \u00d7 24 phase torus (70 mN, Isp 900 s, 150 kg)")
s.shapes.add_picture(os.path.join(A, "phase_torus_dark.png"), Inches(0.4), Inches(1.4), Inches(6.9))
bullets(s, 7.6, 1.6, 5.4, 5.6, [
    (0, [("Vertical stripes: ", GOLD, True), ("t_f is almost constant down a column, so where you meet the tulip dominates", WHITE, False)]),
    (0, [("Range: ", GOLD, True), (f"{float(fast['tf']):.2f} to {float(slow['tf']):.2f} d over the phase torus", WHITE, False)]),
    (0, [("Bottom row and blocks: ", GOLD, True), ("sharp steps are hand-overs between solution families, not noise", WHITE, False)]),
    (0, [("Every cell: ", GOLD, True), ("certified root, second-order sweep clean (0 interior conjugate crossings)", WHITE, False)]),
    (0, [("Use: ", GOLD, True), ("pick the departure and arrival phase and the library returns the costates to fly it", WHITE, False)]),
], size=17)

# ================================================================ CR3BP
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
bullets(s, 0.5, 5.45, 6.4, 1.4, [
    (0, [("DRO: ", GOLD, True), ("stable, retrograde, planar, around the Moon", WHITE, False)]),
    (0, [("Tulip: ", GOLD, True), ("7-petal 3D resonant orbit over the lunar poles", WHITE, False)]),
    (0, [("Orbits from pumpkyn's catalogued families; ", GOLD, True), ("positions set by phase fractions s_D, s_A ∈ [0,1)", WHITE, False)]),
], size=16)
s.shapes.add_picture(os.path.join(A, "cr3bp_geometry.png"), Inches(7.0), Inches(1.45), Inches(6.0))

# ================================================================ min time
s = new_slide(prs, "Minimum-time transfer: thrust all the way, point along the primer",
              "A 70 mN electric thruster moves a 150 kg spacecraft from a lunar DRO to a 7-petal tulip in 16-22 days")
eq1 = eq_png("eq_mintime.png", [
    r"$\min\ t_f \quad \mathrm{s.t.}\quad \dot r = v,\ \ \dot v = g_{\mathrm{CR3BP}}(r,v) + \frac{T}{m}\,\hat u,\ \ \dot m = -\frac{T}{c}$",
    r"$\hat u^{\ast} = -\lambda_v/\|\lambda_v\|,\qquad \|u\| = 1\ \ \mathrm{(all\ burn)}$",
    r"$H(t_f) = 0,\qquad \mathrm{unknowns}\ \ z_8 = [\lambda_r;\ \lambda_v;\ \lambda_m;\ t_f]$",
], fs=19, w=9.6, h=2.1)
s.shapes.add_picture(eq1, Inches(0.5), Inches(1.45), Inches(7.2))
# Mike's comment: make g_CR3BP explicit (gravity + centrifugal + Coriolis)
eqg = eq_png("eq_gcr3bp.png", [
    r"$g_{\mathrm{CR3BP}}(r,v) = \nabla\Omega(r) + 2\,(v_y,\,-v_x,\,0)^{\mathsf{T}}$:",
    r"$\quad g_x = x + 2v_y - (1-\mu)\,\frac{x+\mu}{r_1^3} - \mu\,\frac{x-1+\mu}{r_2^3}$",
    r"$\quad g_y = y - 2v_x - (1-\mu)\,\frac{y}{r_1^3} - \mu\,\frac{y}{r_2^3},\qquad g_z = -(1-\mu)\,\frac{z}{r_1^3} - \mu\,\frac{z}{r_2^3}$",
    r"$\quad r_1 = \|r - (-\mu,0,0)\|\ \mathrm{(Earth)},\qquad r_2 = \|r - (1-\mu,0,0)\|\ \mathrm{(Moon)}$",
], fs=16, w=9.6, h=2.1)
s.shapes.add_picture(eqg, Inches(0.5), Inches(3.05), Inches(7.2))
table(s, 8.1, 1.5, 4.8, [
    ["Engine / spacecraft", "Value"],
    ["Thrust T", "70 mN"],
    ["Specific impulse", "900 s  (c = 8.83 km/s)"],
    ["Initial mass m₀", "150 kg"],
    ["Acceleration T/m₀", "0.47 mm/s²"],
    ["Propellant flow", "0.69 kg/day"],
], colw=[2.3, 2.5], size=15)
bullets(s, 0.5, 4.73, 7.3, 2.6, [
    (0, [("Solved in two stages: ", GOLD, True), ("direct collocation (CasADi/IPOPT) finds the basin; multiple shooting on the PMP field polishes the costates to ~1e-13", WHITE, False)]),
    (0, [("Pontryagin: ", GOLD, True), ("the control is fixed by the costates; the problem becomes an 8-unknown two-point BVP", WHITE, False)]),
    (0, [("Fuel follows time: ", GOLD, True), ("all-burn (a theorem for min time) means propellant = 0.69 kg/day × t_f, so across the library the fastest transfer is also the cheapest", WHITE, False)]),
], size=16)
s.shapes.add_picture(os.path.join(HERE, "..", "indirect", "results", "transfer_3d_anchor.png"),
                     Inches(8.1), Inches(4.3), Inches(4.8))

# ================================================================ APPENDIX (optional)
if APPENDIX:
    # A0 divider
    s = black_slide(prs)
    text(s, 0.9, 2.6, 11.5, 1.0, "Appendix", 44, WHITE, True)
    rule(s, 0.95, 3.6, 2.2)
    text(s, 0.9, 3.8, 11.5, 1.6,
         "The optimality checks behind \"certified\": what each one tests, how it is computed, "
         "and what the whole stack does and does not let us claim", 20, GREY)

    # A1 the logic
    s = new_slide(prs, "Certified = an extremal that also meets a sufficiency theorem",
                  "Problem: fixed departure and arrival states (phases held fixed), free final mass, free t_f, cost t_f",
                  tag="Appendix A1")
    eqA = eq_png("eq_hamiltonian.png", [
        r"$H = 1 + \lambda_r\cdot v + \lambda_v\cdot g_{\mathrm{CR3BP}}(r,v) - s\,T\,Q_{mt},\qquad Q_{mt} = \frac{\|\lambda_v\|}{m} + \frac{\lambda_m}{c}$",
        r"$\mathrm{control:\ throttle}\ s\in[0,1],\ \mathrm{direction}\ \alpha\in S^2;\qquad \dot\lambda = -\partial H/\partial x$",
    ], fs=18, w=10.5, h=1.3)
    s.shapes.add_picture(eqA, Inches(0.5), Inches(1.45), Inches(8.1))
    bullets(s, 0.5, 2.6, 7.9, 4.8, [
        (0, [("Necessary (Pontryagin, first order): ", GOLD, True), ("the root is an extremal: state and costate equations, minimum principle, transversality, H = 0. Any minimizer satisfies them, but so do saddles and maxima", WHITE, False)]),
        (1, [("measured: at a fold of the 70 mN arrival-phase sheet, 12 of 14 candidates passed every first-order check and the conjugate test refuted them", WHITE, False)]),
        (0, [("Sufficient (Bonnard-Caillau-Trélat 2007): ", GOLD, True), ("a normal extremal with the strengthened Legendre condition and no conjugate time in (0, t_f] is a strict strong local minimizer", WHITE, False)]),
        (0, [("Two reductions put us in the theorem's setting: ", GOLD, True), ("strict bang (Q_mt > 0) removes the throttle; then m(t) = 1 − Tt/c exactly, leaving a smooth 6-state problem in α", WHITE, False)]),
        (0, [("Necessary lines are not redundant: ", GOLD, True), ("being an extremal is a hypothesis of the theorem, and the conjugate test is only meaningful along one", WHITE, False)]),
    ], size=15)
    table(s, 8.9, 2.6, 4.1, [
        ["Line state", "Meaning"],
        ["PASS", "checked, within tolerance"],
        ["FAIL", "checked, out of tolerance"],
        ["NOT CHECKED", "never computed: blocks"],
        ["UNRESOLVED", "ran, could not decide: blocks"],
    ], colw=[1.6, 2.5], size=12, rowh=0.38)
    text(s, 8.9, 4.7, 4.1, 2.3,
         "A group passes only if EVERY line is PASS. \"All checks that ran passed\" is vacuously "
         "true when none ran, so an unchecked line blocks the claim. Every certified entry "
         "passes the full stack: certify_root, then a fail-closed audit that re-flies it "
         "from the catalog keys alone.", 13, GREY)

    # A2 necessary
    s = new_slide(prs, "Necessary conditions: six independent first-order checks",
                  "Evaluated on the multiple-shooting solution and on a fresh flight from z₈ alone (certify_root, pmp_pointwise_checks)",
                  tag="Appendix A2")
    table(s, 0.4, 1.4, 12.5, [
        ["Check", "Condition", "How it is computed", "Gate"],
        ["N1  Shooting residual", "costate ODE continuity at K junctions, (r,v)(t_f) = target, λ_m(t_f) = 0, H(t_f) = 0",
         "ms_tfmin damped Newton on the stacked residual R(z)", "‖R‖ ≤ 3e-11"],
        ["N2  Hamiltonian", "H = 1 + λ·f ≡ 0 on [0, t_f] (autonomous, free t_f)",
         "H evaluated at every sample of the flight, not only at t_f", "|H| ≤ 1e-6"],
        ["N3  Flown arrival", "the extremal from z₈ alone reaches the tulip state",
         "fly 16-22 d from z₈, no junction resets; miss in position AND velocity", "≤ 100 km, 10 m/s"],
        ["N4  Transversality", "λ_m(t_f) = 0 (final mass free)",
         "read on the flight; loose vs tight integration gives its uncertainty", "|λ_m(t_f)| ≤ 1e-6"],
        ["N5  Adjoint equations", "λ̇ = −∂H/∂x",
         "central differences of H in the STATE at fixed costate, two steps, Richardson-combined", "rel. err ≤ 1e-7"],
        ["N6  Minimum principle", "applied control minimizes H over s ∈ [0,1], α ∈ S²",
         "control recovered from the field (powered minus coasting; throttle again from the mass row); exact gap below",
         "gap ≤ 1e-12"],
    ], colw=[2.3, 3.9, 4.6, 1.7], size=12, rowh=0.55,
        colors=[WHITE, WHITE, BLUE, GREEN])
    eqN6 = eq_png("eq_gap.png", [
        r"$\mathrm{gap} = H(u,\alpha) - \min_{s,\beta}H = \frac{T}{m}\left(\lambda_v\cdot b + \|b\|\,\|\lambda_v\|\right) + T\left(\max(Q_{mt},0) - \|b\|\,Q_{mt}\right),\quad b = u\,\alpha$",
    ], fs=17, w=12.0, h=0.7)
    s.shapes.add_picture(eqN6, Inches(0.4), Inches(5.55), Inches(10.3))
    text(s, 0.45, 6.4, 12.4, 0.9,
         "Why N5-N6 exist: the shooting residual only says the pieces match each other. A propagator that minimizes the "
         "WRONG Hamiltonian still gives a small residual; the gap, computed from the control actually applied (never from "
         "the formula −λ_v/|λ_v| itself, which would test nothing), catches it. Each check's unit test injects a wrong field and watches it fail.",
         12, GREY)

    # A3 sufficiency
    s = new_slide(prs, "Sufficiency: each hypothesis mapped to a computed margin",
                  "Together with N1-N6 these give a strict strong local minimizer for the fixed phases",
                  tag="Appendix A3")
    table(s, 0.4, 1.4, 12.5, [
        ["Hypothesis", "Mathematical statement", "Our quantity / instrument", "Gate / record"],
        ["S1  Normality", "no abnormal (λ₀ = 0) costate generates the same trajectory",
         "dim S = 1: null space of the lift constraint matrix; Eckart-Young margin σ₆ / measured error at two integration tolerances",
         "margin ≥ 10 (worst 11×)"],
        ["S2  Strict bang", "throttle s = 1 is the unique minimizer with margin: Q_mt > 0",
         "min over the dense flight of Q_mt = |λ_v|/m + λ_m/c", "≥ 1e-5 floor"],
        ["S3  Strengthened Legendre", "H_αα on T_αS² positive definite; here H_αα = (T/m)|λ_v| I₂",
         "min |λ_v| over the dense flight (slope bound |d|λ_v|/dt| ≤ |λ_r| between samples)", "≥ 1e-5 floor"],
        ["S4  No conjugate time", "the Jacobian of the extremal flow in the costate directions keeps full rank on (0, t_f]",
         "sign of det[Φ_rv P, f_rv] at K junctions (ms_conjugate_test) + dense σ₆ scan (conj_spectrum), dips refined 4× and 16×",
         "no sign change; 0 zeros"],
        ["H6  Instrument validity", "λ_m(0) < c/T",
         "h6_margin: (c/T)/λ_m(0), clearance above the Hamiltonian residual", "margin > 1 (worst 5.3×)"],
    ], colw=[2.3, 3.6, 4.8, 1.8], size=12, rowh=0.62,
        colors=[WHITE, WHITE, BLUE, GREEN])
    bullets(s, 0.4, 5.35, 12.5, 2.1, [
        (0, [("S4 in words: ", GOLD, True), ("perturb the initial costate in the 5 directions that keep the problem normalized and watch where the (r,v) state goes; a conjugate time is where some combination of those perturbations (plus a shift along the flow) returns to zero displacement. Past it, a neighbouring extremal reaches the same point and the arc stops being minimizing", WHITE, False)]),
        (0, [("H6 is ours, not the theorem's: ", GOLD, True), ("the reduced Hamiltonian is not conserved, so det = 0 could mean h(t) = 0 instead of a rank drop; since λ_m falls monotonically to 0 and h vanishes only at λ_m = c/T, λ_m(0) < c/T excludes it", WHITE, False)]),
    ], size=13)

    # A4 cross-checks and scope
    s = new_slide(prs, "Cross-checks guard the code; the claim is local and fixed-phase",
                  "Cross-checks are not conditions of the theory, but a failed one still blocks certification",
                  tag="Appendix A4")
    table(s, 0.4, 1.4, 7.4, [
        ["Cross-check", "What it guards against", "Gate"],
        ["X1  pumpkyn tfMin witness", "our solver's bugs: Darin's independent solver must return the stored z₈ unchanged", "|Δz| ≤ 1e-6"],
        ["X2  Independent field", "a wrong equation of motion: pumpkyn's field vs an independently written CR3BP field, state and adjoint rows", "rel. ≤ 1e-10"],
        ["X3  Flight admissibility", "a returned array that is not a flight: reaches t_f, finite, exact all-burn mass law, clear of Moon (1900 km) and Earth (6600 km)", "all hold"],
        ["X4  Fail-closed audit", "packaging errors: every catalog entry re-derived from its own keys and re-flown", "576 OK / 0 bad"],
    ], colw=[2.1, 3.9, 1.4], size=12, rowh=0.75, colors=[WHITE, WHITE, GREEN])
    bullets(s, 8.2, 1.4, 4.9, 5.9, [
        (0, [("What we claim: ", GOLD, True), ("each entry is a strict strong local minimizer among transfers with the same endpoints and the same phases", WHITE, False)]),
        (0, [("Not global: ", GOLD, True), ("another root in the same cell could be faster; the library keeps the fastest one found", WHITE, False)]),
        (0, [("Not phase-optimal: ", GOLD, True), ("phases are held fixed; optimizing them needs its own transversality conditions", WHITE, False)]),
        (0, [("Partly sampled: ", GOLD, True), ("S2, S3 and S4 are tested at finitely many times, so a double zero inside one sub-segment could be missed", WHITE, False)]),
        (0, [("What would close that: ", GOLD, True), ("validated (interval) integration of the variational equations, or an independent Morse-index count", WHITE, False)]),
    ], size=14)

out = os.path.join(HERE, "DRO_tulip_overview_v3.pptx")
prs.save(out)
for k, sl in enumerate(prs.slides, 1):
    print(k, sorted({str(sh.shape_type) for sh in sl.shapes}))
print("saved", out)
