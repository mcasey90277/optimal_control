# sundman_minfuel/attic — superseded drivers (do not use)

Atticked in cleanup Phase 1 (2026-07-09). All replaced by the canonical
toolchain `minfuel_config` / `minfuel_at_tf` / `aggregate_front` +
`orchestrate/*.sh`; they also reference pre-migration result filenames
(energy_<%.2f>.mat etc.) that no longer exist, so they are broken as-is.
`tf_front.png` / `tf_dv_front.png` are the historical scatter plots
(non-extremal local minima — see HONEST_EVALUATION_DV_TF_FRONT.md).
