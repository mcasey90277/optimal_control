# HALO_HALO — TODO

- [ ] **Full L1↔L2 catalog** (2 L1 × 11 L2 members × phasing × rungs).
  Blocked on the arrival-period schema axis (see README); implement the
  SDD's "new parameter axis" registration in
  `costate_common/build_costate_catalog_family` + the pickers, then this
  campaign is a front-door parameterization like halo/DPO.
- [ ] **Both directions:** L2 → L1 is a distinct problem (no symmetry
  argument covers it); probe one pair before assuming the ladder behaves.
- [ ] **Phasing:** probes ran a single (0, 0) phase cell; the catalog needs
  the 6×6 torus.
- [ ] **10 N NaN watch item** (far pair): tfMin acceptance returned NaN
  with a healthy 0.08 km flown miss — diagnose if it recurs at catalog
  scale (likely a tfMin exception, possibly near a family wall).
- [ ] Northern-branch mirror statement (pm = +1) once a deliverable exists.
