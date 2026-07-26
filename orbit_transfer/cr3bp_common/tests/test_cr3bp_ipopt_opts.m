% TEST_CR3BP_IPOPT_OPTS  Byte-identity gate for the extracted IPOPT option set.
%
% The extraction of cr3bp_ipopt_opts.m removed ~20 duplicated option lines from
% each of three CERTIFIED solvers. The refactor is only safe if the helper
% reproduces what each of them built inline, EXACTLY -- a different mu_init or
% a missing warm_start push would change convergence, silently, on solvers
% whose results are certified.
%
% So this test does not check that the helper "looks right". It reconstructs
% the pre-extraction inline struct verbatim (transcribed from the solvers as
% they stood at commit ece0680) and asserts isequal against the helper's
% output, for BOTH warm-start regimes. It is the gate the extraction rests on,
% and it is why no re-solve of a certified row was needed to justify the change.
%
% If this test ever fails, the helper and the recorded historical behaviour
% have diverged: decide deliberately which is correct, and if the helper is,
% update the reference below WITH a note on why the behaviour changed.
%
% INPUTS:  none
% OUTPUTS: none (throws on mismatch)
%
% REFERENCES:
%   [1] cr3bp_common/cr3bp_ipopt_opts.m (the helper under test).
%   [2] orbit_transfer/CODE_STRUCTURE.md (Tier-0 plan and the gate rule).
here = fileparts(mfilename('fullpath'));
addpath(fullfile(here, '..'));

for tf = [true false]
    maxIter = 3000;

    % --- reference: the inline struct as the three solvers built it ----------
    ref = struct();
    ref.print_time      = true;
    ref.ipopt.max_iter  = maxIter;
    ref.ipopt.tol       = 1e-7;
    ref.ipopt.nlp_scaling_method = 'gradient-based';
    ref.ipopt.linear_solver      = 'mumps';
    ref.ipopt.print_level        = 5;
    if tf
        ref.ipopt.mu_strategy                 = 'monotone';
        ref.ipopt.warm_start_init_point       = 'yes';
        ref.ipopt.mu_init                     = 1e-4;
        ref.ipopt.warm_start_bound_push       = 1e-9;
        ref.ipopt.warm_start_bound_frac       = 1e-9;
        ref.ipopt.warm_start_slack_bound_push = 1e-9;
        ref.ipopt.warm_start_slack_bound_frac = 1e-9;
        ref.ipopt.warm_start_mult_bound_push  = 1e-9;
    else
        ref.ipopt.mu_strategy           = 'monotone';
        ref.ipopt.warm_start_init_point = 'yes';
        ref.ipopt.mu_init               = 0.1;
    end

    got = cr3bp_ipopt_opts(maxIter, tf);

    % Field sets must match exactly -- a MISSING option is the dangerous case
    % (IPOPT would silently use its own default), and an EXTRA one changes
    % behaviour just as much.
    fr = sort(fieldnames(ref.ipopt));  fg = sort(fieldnames(got.ipopt));
    assert(isequal(fr, fg), ...
        ['test_cr3bp_ipopt_opts (warmTight=%d): ipopt field sets differ.\n' ...
         '  only in reference: %s\n  only in helper   : %s'], tf, ...
        strjoin(setdiff(fr, fg), ', '), strjoin(setdiff(fg, fr), ', '));

    for k = 1:numel(fr)
        f = fr{k};
        assert(isequal(ref.ipopt.(f), got.ipopt.(f)), ...
            'test_cr3bp_ipopt_opts (warmTight=%d): ipopt.%s differs', tf, f);
    end
    assert(isequal(ref.print_time, got.print_time), 'print_time differs');
    assert(isequal(ref, got), ...
        'test_cr3bp_ipopt_opts (warmTight=%d): structs differ beyond .ipopt', tf);

    fprintf('  warmTight=%d : %2d ipopt options, all identical to the pre-extraction inline struct\n', ...
            tf, numel(fr));
end

% maxIter must be threaded, not hardcoded -- a stuck cap would silently
% under-iterate the deep rungs, which is a documented failure mode in this repo.
p1 = cr3bp_ipopt_opts(123, true);
p2 = cr3bp_ipopt_opts(456, true);
assert(p1.ipopt.max_iter == 123 && p2.ipopt.max_iter == 456, ...
    'test_cr3bp_ipopt_opts: maxIter is not threaded through');

fprintf('test_cr3bp_ipopt_opts: ALL PASS\n');
