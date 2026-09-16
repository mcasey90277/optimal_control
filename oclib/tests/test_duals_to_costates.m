function ok = test_duals_to_costates()
%% Purpose:
%
%   Unit test for oc.duals_to_costates on synthetic multipliers, so it runs
%   from a fresh clone (the function's self-demo needs a gitignored campaign
%   file). Checks the contract the header states, one rule at a time:
%     1. STATIONS: Hermite-Simpson -> interval midpoints, trapezoid -> left
%        nodes, trapezoid-nodal -> all N+1 nodes, on a NON-uniform grid.
%     2. TRAPEZOID-NODAL AVERAGE: interior columns are the step-weighted
%        adjacent-interval average, endpoints one-sided.
%     3. DEFECT SCALE: 'h' divides each column by its interval width.
%     4. SIGN VOTE: multipliers whose primer -lambda_v aligns with the
%        control keep sign +1; the flipped multipliers are put back (-1).
%     5. LAMBDA_T: a lifted-t_f row resolves to +1 after the sign.
%     6. NO CONTROL: no uDir -> sign left as-is, NaN margin, a note.
%     7. ACTIVE FLAGS: the active fraction is reported.
%     8. UNIFORM MESH from tf when tNodes is empty.
%     9. REFUSALS: unknown scheme, tNodes/mu size mismatch, no tf.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 09/16/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));                         % oclib root -> +oc visible
ok = true;
rngWas = rng(7);                                  % reproducible fixture
restoreRng = onCleanup(@() rng(rngWas));

%% Synthetic multipliers on a non-uniform grid:
M     = 12;
tN    = cumsum([0, 0.5 + rand(1, M)]);            % strictly increasing
h     = diff(tN);
uDir  = randn(3, M);
uDir  = uDir ./ sqrt(sum(uDir.^2, 1));            % unit control directions
lamV  = -2.5 * uDir;                              % primer -lam_v || u
mu    = [randn(3, M); lamV; randn(1, M)];         % [lam_r; lam_v; lam_m]
lamTf = ones(1, M) + 1e-6*randn(1, M);

%% 1. Stations per scheme:
[~, tS] = oc.duals_to_costates(struct('scheme', 'hermite-simpson', 'mu', mu, 'tNodes', tN));
ok = check('hermite-simpson stations are interval midpoints', ...
           max(abs(tS - (tN(1:end-1) + h/2))) < 1e-14) && ok;
[lamT, tS] = oc.duals_to_costates(struct('scheme', 'trapezoid', 'mu', mu, 'tNodes', tN));
ok = check('trapezoid stations are left nodes, mu unchanged', ...
           isequal(tS, tN(1:end-1)) && isequal(lamT, mu)) && ok;
[lamN, tS] = oc.duals_to_costates(struct('scheme', 'trapezoid-nodal', 'mu', mu, 'tNodes', tN));
ok = check('trapezoid-nodal stations are all N+1 nodes', ...
           isequal(tS, tN) && isequal(size(lamN), [7 M+1])) && ok;

%% 2. Trapezoid-nodal weighted average:
kn  = 5;
ref = (h(kn-1)*mu(:,kn-1) + h(kn)*mu(:,kn)) / (h(kn-1) + h(kn));
ok = check('interior node is the step-weighted adjacent average', ...
           max(abs(lamN(:,kn) - ref)) < 1e-14) && ok;
ok = check('endpoint columns are one-sided', ...
           isequal(lamN(:,1), mu(:,1)) && isequal(lamN(:,M+1), mu(:,M))) && ok;

%% 3. Defect scale 'h':
lamH = oc.duals_to_costates(struct('scheme', 'trapezoid', 'mu', mu, 'tNodes', tN, ...
                                   'defectScale', 'h'));
ok = check('''h'' scaling divides each column by its interval width', ...
           max(max(abs(lamH - mu./h))) < 1e-14) && ok;

%% 4. Sign vote, both orientations:
s = struct('scheme', 'hermite-simpson', 'mu', mu, 'tNodes', tN, 'uDir', uDir, ...
           'lamTf', lamTf);
[lamP, ~, dP] = oc.duals_to_costates(s);
ok = check(sprintf('aligned multipliers keep sign +1 (margin %.2f)', dP.voteMargin), ...
           dP.sign == 1 && dP.voteMargin == 1 && isequal(lamP, mu)) && ok;
s.mu = -mu;  s.lamTf = -lamTf;
[lamF, ~, dF] = oc.duals_to_costates(s);
ok = check('flipped multipliers are put back (sign -1, lam == original)', ...
           dF.sign == -1 && dF.voteMargin == 1 && isequal(lamF, mu)) && ok;

%% 5. Lambda_t after the sign:
ok = check(sprintf('lambda_t resolves to +1 after the flip (%.6f)', dF.lamT), ...
           dF.lamTOK && abs(dF.lamT - 1) < 1e-3) && ok;
s.lamTf = -2*lamTf;
[~, ~, dBad] = oc.duals_to_costates(s);
ok = check('a lambda_t far from +1 is flagged, with a note', ...
           ~dBad.lamTOK && any(contains(dBad.notes, 'lambda_t'))) && ok;

%% 6. No control supplied:
[lam0, ~, d0] = oc.duals_to_costates(struct('scheme', 'trapezoid', 'mu', -mu, 'tNodes', tN));
ok = check('no uDir: sign unresolved, NaN margin, a note says so', ...
           d0.sign == 1 && isnan(d0.voteMargin) && isequal(lam0, -mu) ...
           && any(contains(d0.notes, 'sign NOT resolved'))) && ok;

%% 7. Active path-constraint flags:
af = false(1, M);  af(1:3) = true;
[~, ~, dA] = oc.duals_to_costates(struct('scheme', 'trapezoid', 'mu', mu, 'tNodes', tN, ...
                                         'activeFlags', af));
ok = check('active fraction reported (3 of 12)', dA.activeFraction == 0.25) && ok;

%% 8. Uniform mesh from tf:
[~, tU] = oc.duals_to_costates(struct('scheme', 'trapezoid', 'mu', mu, 'tNodes', [], 'tf', 6));
tRef = linspace(0, 6, M+1);
ok = check('empty tNodes + tf builds the uniform mesh', ...
           isequal(tU, tRef(1:M))) && ok;

%% 9. Refusals:
ok = check('unknown scheme refused', throwsId(@() oc.duals_to_costates( ...
           struct('scheme', 'euler', 'mu', mu, 'tNodes', tN)), 'duals_to_costates:scheme')) && ok;
ok = check('tNodes/mu size mismatch refused', throwsId(@() oc.duals_to_costates( ...
           struct('scheme', 'trapezoid', 'mu', mu, 'tNodes', tN(1:end-1))), 'duals_to_costates:size')) && ok;
ok = check('empty tNodes without tf refused', throwsId(@() oc.duals_to_costates( ...
           struct('scheme', 'trapezoid', 'mu', mu, 'tNodes', [])), 'duals_to_costates:tf')) && ok;

fprintf('TEST_DUALS_TO_COSTATES: %s\n', passText(ok));
end

% ------------------------------------------------------------------------
function ok = check(name, cond)
%% Purpose:
%
%   Print one gate line and return its verdict.
%
ok = logical(cond);
if ok, tag = 'PASS'; else, tag = 'FAIL'; end
fprintf('  [%s] %s\n', tag, name);
end

function tf = throwsId(fn, id)
%% Purpose:
%
%   True when fn() throws an error whose identifier is id.
%
tf = false;
try
    fn();
catch ME
    tf = strcmp(ME.identifier, id);
end
end

function s = passText(ok)
%% Purpose:
%
%   'ALL PASS' or 'FAILED'.
%
if ok, s = 'ALL PASS'; else, s = 'FAILED'; end
end
