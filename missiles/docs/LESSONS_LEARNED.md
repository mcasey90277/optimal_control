# Lessons Learned — Missile Trajectory Library

Running log. Newest entries at the top. Record what broke, what fixed it, and
what a future reader would otherwise rediscover the hard way.

## 2026-08-06 — Vacuum energy conservation is the test that matters

The glide EOM has six equations and many chances for a sign error. Propagating
a vacuum ballistic arc and checking that V^2/2 - mu/r holds constant to 1e-8
catches nearly all of them in one assertion, because any sign flip in gravity
or the centrifugal term breaks conservation immediately. Approximate checks
against analytic glide solutions do not — they tolerate a percent of error by
construction, which is enough to hide a wrong term.

## 2026-08-06 — `-batch` needs `run(...)`, not a bare relative path

Invoking the test harness as `matlab -batch "cd(...); tests/run_tests"` does
NOT work: MATLAB parses `tests/run_tests` as the expression `tests / run_tests`
(division of two undefined names) and errors with `Unrecognized function or
variable 'tests'` before `run_tests.m` ever executes. The working forms are
`matlab -batch "cd(...); run('tests/run_tests')"` or
`matlab -batch "cd(...); addpath('tests'); run_tests"`. Use `run('tests/run_tests')`
for all future headless invocations of this harness.

## 2026-08-06 — Constants live in one place

`missileConst` is the single source of truth for Earth and air constants.
pumpkyn's `getConst` was not used for these: it carries `g` and `R` but no
`muE`, no Earth radius, and two fields known to be wrong (`deg2ArcSec`, `c`).
Mixing the two would make it ambiguous which constant a routine actually used.
