### Task 12: Theory note (`doc/booster_landing_note.tex`)

**Files:**
- Create: `doc/booster_landing_note.tex`

**Interfaces:**
- Consumes: `results/` figures + the flagship numbers (real, from the runs — never invented).
- Produces: compiled PDF. Sections:
  1. **Problem** — 3-DOF PDG statement, F9-class parameters table, the hoverslam T/W argument
  2. **Two routes** — HS collocation NLP; lossless convexification (change of variables, the relaxation, statement of the Açıkmeşe–Ploen tightness theorem with proof sketch via the maximum principle; Taylor mass-bound construction)
  3. **Certification** — the five gates, table of flagship values (from `booster_run.mat`)
  4. **Tracking** — TVLQR derivation sketch, saturation-near-terminal-arc discussion (expected physics, shown honestly)
  5. **Monte Carlo** — dispersions table, footprint figure, success rate
  6. **Phase 2** — drag model, Δfuel result, what drag-free guidance misses
  7. **Future work** — 6-DOF, full return profile, SCvx, real conic solver
- References: Açıkmeşe & Ploen 2007; Blackmore, Açıkmeşe & Scharf 2010; Kelly 2017; Anderson & Moore. All real — no invented citations.

- [ ] **Step 1: Write the note**

Standard article class, house LaTeX conventions. Pull every quoted number from `results/booster_run.mat` (load it, read the fields) — if a number in the draft can't be traced to the .mat or a test output, delete it.

- [ ] **Step 2: Compile + clean aux**

Run: `cd /Users/msc/Desktop/optimal_control/booster_landing/doc && /Library/TeX/texbin/pdflatex -interaction=nonstopmode booster_landing_note.tex && /Library/TeX/texbin/pdflatex -interaction=nonstopmode booster_landing_note.tex && rm -f *.aux *.log *.out *.toc *.lof *.lot *.fls *.fdb_latexmk *.synctex.gz`
Expected: clean compile, no undefined references.

- [ ] **Step 3: Verify figures/citations**

Run: `~/ai_council/venv/bin/python ~/Documents/myLatex/tools/paper_verify.py /Users/msc/Desktop/optimal_control/booster_landing/doc/booster_landing_note.tex --deep-figs`
Expected: green — no FAIL on cites/figures (house rule: mandatory before calling doc work done).

- [ ] **Step 4: Commit**

```bash
cd /Users/msc/Desktop/optimal_control
git add booster_landing/doc/booster_landing_note.tex
git commit -m "booster_landing: theory note (PDG, convexification, gates, TVLQR, MC)

Co-Authored-By: Claude Fable 5 <noreply@anthropic.com>"
```

(Then offer Mike the doc-review pass on the note — decided at spec time to review the note, not the spec.)

---

## Execution order & checkpoints

Tasks are sequential (each consumes the previous interfaces). Natural review checkpoints for Mike: after Task 5 (the scientific core — two methods, one answer), after Task 10 (flagship), after Task 12 (note).

## Self-review notes (done at write time)

- Spec coverage: every spec section maps to a task (problem/params→1–2, solvers→3–4, certify→5, tracking→6–7, MC→8, viz→9, front door→10, Phase 2→11, note→12, tests woven throughout, pumpkyn conventions in Global Constraints).
- Type consistency: `sol.{t,tf,mf,X,U}` identical meaning across both solvers (convex converts back to thrust-in-N); `ctrl.{tgrid,K,xnom,Tnom}` consumed by sim exactly as produced; `P` merge pattern uniform.
- Known judgment calls flagged inline: G5 primer sign (flip once, don't abs), coarse-grid G3 tolerance scaling, golden-section unimodality assumption, drag-case primer law.
