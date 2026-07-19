function AS = sosc_active_set(R, K, tol)
% SOSC_ACTIVE_SET  Classify inequality rows active/strong/weak, assemble the
% active Jacobian (equalities + strongly-active inequalities), flag LICQ.
%
% INPUTS: R - recover struct; K - kkt_residual struct (for sign); tol - defaults
% OUTPUTS: AS - struct .A[sparse m_active x n] .m_active .nEq .nStrong .nWeak
%               .weakLabels{cell} .licq
% REFERENCES: process/DESIGN_sosc.md sec 4.4.
muThresh = tol.mu * max(1, max(abs(R.lam_g)));
eqRows = []; strongRows = []; nWeak = 0; weakLabels = {};
for i = 1:numel(R.creg)
    c = R.creg(i);
    if strcmp(c.kind,'eq'), eqRows = [eqRows, c.rows]; continue; end %#ok<AGROW>
    for j = 1:numel(c.rows)
        r = c.rows(j);
        gv = R.gval(r);
        if strcmp(c.kind,'ineqHi'), slack = c.bound - gv; else, slack = gv - c.bound; end
        if slack < tol.active                       % active
            if abs(R.lam_g(r)) > muThresh
                strongRows = [strongRows, r]; %#ok<AGROW>
            else
                nWeak = nWeak + 1;
                nd = ''; if ~isempty(c.node), nd = sprintf(', node %d', c.node(min(j,numel(c.node)))); end
                weakLabels{end+1} = sprintf('%s%s (slack %.1e, lam %.1e)', c.label, nd, slack, R.lam_g(r)); %#ok<AGROW>
            end
        end
    end
end
actRows = sort([eqRows, strongRows]);
AS.A = R.A_all(actRows, :);
AS.m_active = numel(actRows);
AS.nEq = numel(eqRows); AS.nStrong = numel(strongRows); AS.nWeak = nWeak;
AS.weakLabels = weakLabels;
AS.licq = (sprank(AS.A) == AS.m_active);
end
