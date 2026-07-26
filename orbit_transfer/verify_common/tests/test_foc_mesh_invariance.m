% TEST_FOC_MESH_INVARIANCE  Regression test for finding I1 (external review
% 2026-07-25): the regular-switching statistic must be MESH-INVARIANT for a
% physically regular (transversal) switch.
%
% Builds the same bang-bang toy OCP (double integrator, min int u dt, one
% switch) at two mesh densities and asserts that rep.sdotMinRel agrees to
% within a factor of ~2, while the LEGACY un-normalized statistic does not.
% The legacy assertion is the point: it documents the defect this fix removes,
% and will alert us if someone silently reverts the normalization.
%
% REFERENCES:
%   [1] verify_common/foc_check.m section (7) (the corrected statistic).
%   [2] doc/review_2026-07-25_{gpt56terra,gemini31pro}.md (both reviewers
%       derived the same deweight -> divided-difference -> nondimensionalize fix).
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_verify_common;
addpath(fullfile(getenv('HOME'), 'casadi-3.7.0'));  import casadi.*

Ns = [60, 240];                      % 4x refinement
R  = nan(size(Ns));  Rleg = nan(size(Ns));  nsw = nan(size(Ns));
for q = 1:numel(Ns)
    N = Ns(q);  T = 2;  sg = linspace(0,1,N+1).';  ds = diff(sg).';
    opti = casadi.Opti();
    X = opti.variable(2, N+1);  U = opti.variable(1, N+1);
    creg = struct('label',{},'rows',{});
    f = [X(2,:); U(1,:)];
    r0 = size(opti.g,1)+1;
    for k = 1:N
        opti.subject_to(X(:,k+1) - X(:,k) - (ds(k)*T/2)*(f(:,k)+f(:,k+1)) == 0);
    end
    creg(end+1) = struct('label','defect','rows',r0:size(opti.g,1)); %#ok<AGROW>
    r0 = size(opti.g,1)+1;  opti.subject_to(U(1,:).' >= 0);
    creg(end+1) = struct('label','thrLo','rows',r0:size(opti.g,1)); %#ok<AGROW>
    r0 = size(opti.g,1)+1;  opti.subject_to(U(1,:).' <= 1);
    creg(end+1) = struct('label','thrHi','rows',r0:size(opti.g,1)); %#ok<AGROW>
    opti.subject_to(X(:,1) == [0;0]);  opti.subject_to(X(1,end) == 1);
    opti.minimize(sum((ds*T/2).*(U(1,1:N)+U(1,2:N+1))));
    opti.set_initial(X, [linspace(0,1,N+1); 0.5*ones(1,N+1)]);
    opti.set_initial(U, 0.5*ones(1,N+1));
    opti.solver('ipopt', struct('print_time',false), struct('print_level',0));
    sol = opti.solve();
    out = struct('X', full(sol.value(X)), 'U', full(sol.value(U)), ...
                 'model', struct('opti', opti, 'creg', creg));
    rep = foc_check(out, sg, foc_manifest('toy'), struct());
    R(q)   = rep.sdotMinRel;
    nsw(q) = rep.nSwitches;
    % The un-normalised LEGACY statistic is reconstructed HERE rather than
    % carried as a vestigial production field (removed 2026-07-26 once the
    % confound was confirmed). Keeping it in the test preserves the
    % two-directional assertion below, which is what makes this a regression
    % test rather than a tautology.
    burnQ = rep.Sd < 0;
    swQ   = find(diff(burnQ) ~= 0);
    Rleg(q) = min(abs(rep.Sd(swQ+1) - rep.Sd(swQ))) / ...
              max(median(abs(rep.Sd(~burnQ))), 1e-30);
    fprintf('  N=%4d : nSwitch=%d  sdotMinRel=%.4e   (legacy raw %.4e)\n', ...
            Ns(q), nsw(q), R(q), Rleg(q));
end

assert(all(nsw == 1), 'Test: both meshes must find exactly 1 switch (got %s)', mat2str(nsw));

ratio    = max(R)/min(R);
ratioLeg = max(Rleg)/min(Rleg);
assert(ratio < 2.5, ...
    ['Test: normalized sdotMinRel must be mesh-invariant for a regular switch; ' ...
     'got %.4e vs %.4e (ratio %.2f) across a 4x refinement -- the I1 ' ...
     'normalization is broken or reverted'], R(1), R(2), ratio);
assert(ratioLeg > 2.5, ...
    ['Test: the LEGACY statistic was supposed to be mesh-DEPENDENT (that is ' ...
     'the defect I1 fixes); it came back invariant (ratio %.2f), so this test ' ...
     'no longer discriminates and needs rethinking'], ratioLeg);
fprintf(['test_foc_mesh_invariance: ALL PASS -- normalized ratio %.2f (<2.5), ' ...
         'legacy ratio %.1f (mesh-dependent, as expected)\n'], ratio, ratioLeg);
