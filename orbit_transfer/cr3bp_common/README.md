# cr3bp_common — shared CR3BP GTO-transfer library

Single source of truth for the problem definition shared by the `GTO_tulip`
and `GTO_ELFO` campaigns (both Earth–Moon CR3BP, same 15 kg / 25 mN / Isp
2100 s spacecraft). Extracted 2026-07-21 so the campaigns stopped reaching
into each other's folders for these files.

| file | what |
|---|---|
| `cr3bp_lt_params.m` | CR3BP low-thrust dynamics parameters (canonical rotating-frame units, μ*, thrust/Isp nondimensionalization). |
| `minfuel_config.m` | Campaign configuration: thrust/mass/Isp, homotopy schedules, directory map, and BOTH min-time anchors — `cfg.tfMin` (tulip, 6.2906939607 ND) and `cfg.tfMin_elfo` (ELFO, 6.0961534862 ND; factors are labeled against the target's own anchor since the 2026-07-21 rebase). |
| `gto_tulip_endpoints.m` | Tulip problem endpoints (used by tulip direct+indirect AND by elfo). |
| `gto_elfo_endpoints.m` | ELFO problem endpoints. |
| `setup_cr3bp_common.m` | Adds this folder + the pumpkyn toolbox to the path (asserts pumpkyn exists). Called by every GTO module's `setup_paths.m`. |

Deliberately NOT here: `GTO_tulip/direct/PSR/lib/` vendors its own frozen
copies (self-contained by design — see the PSR notes); the earth-GEO
campaigns are 2-body and have their own params.
