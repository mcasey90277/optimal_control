Finite-difference checks passed for `lt_dynamics.m` A/B and `nlp_constraints.m` gradient layout; no CR3BP/Jacobian sign or transpose error found.

- **[EFFICIENCY]** `nlp_constraints.m:50` — `lt_dynamics` is always called with three outputs, so A/B are built even for value-only calls like warm-start defect checks. Call `[F,A,B]` only inside `if nargout > 2`; otherwise call `F = lt_dynamics(...)`.
- **[ROBUSTNESS]** `build_guess.m:95` — Mesh normalization assumes `tauArc(1)==0`, finite values, and strictly increasing samples. Normalize with `(tauArc-tauArc(1))/(tauArc(end)-tauArc(1))`, remove duplicates in `density_matched_mesh`, and reject/sort nonmonotone abscissae before `interp1`.
- **[ROBUSTNESS]** `build_guess.m:66` — Primer/velocity normalization can divide by zero if `lamV` or `Vg` has a near-zero column. Guard norms with a tolerance and use a neighboring/fallback unit direction.
- **[ROBUSTNESS]** `solve_tfmin_nlp.m:31` — Solver trusts `sigma`/`Z0` shape and monotonicity; malformed inputs silently create wrong segment widths or reshape failures. Add upfront assertions for `numel(Z0)==10*(N+1)+1`, `sigma(1)==0`, `sigma(end)==1`, and `all(diff(sigma)>0)`.
- **[STYLE]** `NLP_lowThrust_GTO_Tulip.m:64` — `rvInt` is propagated for one GTO period but only `rvInt(1,:)` is used, making `P0`/propagation dead or indexing-intent ambiguous. If departure is initial GTO, use `rv0 = rv0ND`; if after a coast, use `rvInt(end,:)` or an explicit `idx_0`.
- **[STYLE]** `NLP_lowThrust_GTO_Tulip.m:24` — Header documents `.errPos_km`, `.errVel_kms`, and `.U`, but output uses `.devPos_km`, `.devVel_kms`, and `.W`. Rename fields or fix the header.
- **[STYLE]** `NLP_lowThrust_GTO_Tulip.m:84` — `Z` is assigned but unused. Replace with `[~, nlp] = ...` or include `Z` in `out`.
- **[ROBUSTNESS]** `setup_paths.m:14` — Hard-coded `~/Desktop/proj7/external/pumpkyn/src` makes the solver non-portable. Resolve pumpkyn relative to the project, accept an input path, or document/set it through configuration.

Top 5 priorities:
1. Harden `density_matched_mesh` against bad/duplicate/nonzero-start abscissae.
2. Add `sigma`/`Z0` validation before `fmincon`.
3. Clarify/remove the GTO propagation indexing ambiguity.
4. Avoid unnecessary A/B construction in value-only constraint calls.
5. Guard all unit-vector normalizations against zero norms.
