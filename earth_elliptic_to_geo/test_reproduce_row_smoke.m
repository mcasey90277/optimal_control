% test_reproduce_row_smoke.m — 10 N re-solved FROM SCRATCH (coldB) via the
% keep-best-mass multi-start. Under the one-sided-verify decision (memory
% tenN-minfuel-razor-basin) the gate is MASS: the reproduced solution must be
% at least as good as the campaign floor (and is expected to BEAT it, ~1378.46
% vs 1377.10 kg). Structure (switches/revs) is reported, not asserted.
addpath(fullfile(getenv('HOME'),'casadi-3.7.0'));
here = fileparts(mfilename('fullpath')); cd(here);
row  = reproduce_row(10);          % coldB: no prev rung; writes results/repro/
cert = table3_certified(10);
assert(row.certified == true, '10 N reproduced solution certified');
assert(row.m_f_kg >= cert.m_f_kg - 0.5, ...
    sprintf('10 N reproduced m_f=%.4f must be >= campaign floor %.4f', row.m_f_kg, cert.m_f_kg - 0.5));
assert(isfile(fullfile(module_root(),'results','repro','REPRO_row_T100.mat')), 'REPRO row written');
if row.m_f_kg > cert.m_f_kg + 1e-3
    fprintf('  reproduce_row(10) BEAT the campaign: %.4f vs %.4f kg (+%.4f) [%d sw / %.3f rev]\n', ...
        row.m_f_kg, cert.m_f_kg, row.m_f_kg - cert.m_f_kg, row.switches, row.revs);
end
fprintf('test_reproduce_row_smoke PASSED\n');
