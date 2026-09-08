# PLQ / Huber smoothing — literature check (2026-09-07)

Collected while answering "is our Huber/PLQ energy->fuel homotopy novel?"
(see `orbit_transfer/DRO_tulip/FINDINGS.md` §22–29 and
`orbit_transfer/doc/algorithms_orbit_transfer.tex` §4.2). Verdict of the first
pass: no published indirect low-thrust smoothing family is Huber / piecewise
linear-quadratic, and none keeps a DISCONTINUOUS regularized throttle law
propagated with saltation matrices — but the mathematical objects are known
(eps = elastic net; Huber = Moreau envelope of |s|; huberc = one instance of
the Sandia "any monotone saturation" recipe). Not a Scholar-level review.

| file | what it is / why it is here |
|---|---|
| `Tafazzol_Taheri_2024_AAS24-302_...` | The 2024 comparison of regularizations for indirect min-fuel low-thrust (tanh vs L2-norm). Its reference list is the place to start a proper novelty check. |
| `Wang_etal_2023_new_smoothing_L2norm_bangbang.pdf` | Normalized L2-norm smoothing of the sign of the switching function (arXiv 2309.03069). |
| `Heidrich_Sparapany_Grant_2020_SAND2020-13728C_...` | Sandia: general recipe — pick ANY monotone saturation u~(u) of the switching function, derive the matching error term L~ = ∫ ξ u~'(ξ) dξ. Our huberc two-slope ramp is one instance; a discontinuous law (Huber's jump) is outside its assumptions. |
| `Schneider_Wachsmuth_2018_ESAIM-COCV_...` | L2-regularized L1 (group-sparsity) control cost for ODEs, motivated by a satellite with isotropic thrust — the elastic-net (= our eps family) analysed rigorously: sparsity, bang-bang at alpha = 0, stability, discretization. |
| `Kong_etal_2023_saltation_matrices_hybrid_systems.pdf` | Saltation matrices tutorial (arXiv 2306.06862) — the jump in the STM at a control discontinuity that our Huber propagator needed (FINDINGS §22). |

Paywalled / not fetched (cite from memory or library): Silva & Trélat 2010
(IEEE TAC 55(11) 2488–2499, smooth regularization of bang-bang OCPs — the
convergence theory behind the homotopy line; author copy at
ljll.fr/trelat/fichiers/SilvaTrelat_TAC2010.pdf sits behind a bot wall, fetch
it from a browser); Bertrand & Epenoy 2002
(OCAM 23, quadratic smoothing, the eps family); Taheri & Junkins 2018 (JGCD,
hyperbolic-tangent smoothing) and AAS 19-496 (tanh + STM); the JAS 2023
comparison (Springer, 10.1007/s40295-023-00417-4); Mall & Taheri UTM (JSR).
