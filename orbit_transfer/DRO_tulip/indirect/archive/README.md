# archive — historical one-off scripts

These are **not production code and are not on the supported path.** They are
kept only because they are the record of how a result was first obtained, and
some of the results they produced live in `results/*.mat`, which is gitignored.

They have hard-coded scratchpad paths that no longer resolve. Do not run them;
read them.

| file | what it was | superseded by |
|---|---|---|
| `direct70_2026-09-08.m` | first 70 mN DRO → tulip solve, by direct collocation seeded with the FULL 75.5 mN trajectory (state, control and an exactly consistent mass history, not a bare `lambda0`) | `run_dro_tulip` |
| `certify70_2026-09-08.m` | the covector harvest → `ms_bvp` → `tfMin` acceptance run that turned that solve into the certified anchor, sweeping K = 24/48/96 | `certify_root` / `certify_crossing` |

The supported way to obtain one 70 mN transfer is:

```matlab
T = run_dro_tulip(sD, sA);
```

which serves it from `dro_tulip_library` when it is already solved and
otherwise walks there by continuation, and which puts every answer through
the same gate stack.
