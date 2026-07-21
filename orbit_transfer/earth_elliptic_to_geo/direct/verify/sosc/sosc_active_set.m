function AS = sosc_active_set(R, K, tol)
% SOSC_ACTIVE_SET  Classify inequality rows active/strong/weak, assemble the
% active Jacobian (equalities + strongly-active inequalities), flag LICQ.
% Active/kind are sourced from the canonical NLP bounds R.lbg/R.ubg (§11.2):
% a row is an inequality iff lbg~=ubg, and active iff its distance to the
% nearest finite bound < tol.active. Equalities are always active. Human weak
% labels still come from creg.label/creg.node.
%
% INPUTS: R - recover struct (needs .lbg .ubg); K - kkt_residual struct;
%         tol - defaults
% OUTPUTS: AS - struct .A[sparse m_active x n] .m_active .nEq .nStrong .nWeak
%               .weakLabels{cell} .licq
% REFERENCES: process/DESIGN_sosc.md sec 4.4, 11.2.
lbg = R.lbg(:);  ubg = R.ubg(:);
isEq = (lbg == ubg);

% distance to the nearest finite bound (inf where a side is unbounded)
distLo = R.gval - lbg;  distLo(~isfinite(lbg)) = inf;
distHi = ubg - R.gval;  distHi(~isfinite(ubg)) = inf;
slackRow = min(distLo, distHi);

% per-row human label/node from creg (for weak-node reporting only)
rowLabel = strings(R.m,1);  rowNode = nan(R.m,1);
for i = 1:numel(R.creg)
    c = R.creg(i);
    for j = 1:numel(c.rows)
        r = c.rows(j);
        rowLabel(r) = c.label;
        if ~isempty(c.node), rowNode(r) = c.node(min(j,numel(c.node))); end
    end
end

muThresh = tol.mu * max(1, max(abs(R.lam_g)));
eqRows = []; strongRows = []; nWeak = 0; weakLabels = {};
for r = 1:R.m
    if isEq(r), eqRows(end+1) = r; continue; end %#ok<AGROW>
    if slackRow(r) < tol.active                     % active inequality
        if abs(R.lam_g(r)) > muThresh
            strongRows(end+1) = r; %#ok<AGROW>
        else
            nWeak = nWeak + 1;
            nd = ''; if ~isnan(rowNode(r)), nd = sprintf(', node %d', rowNode(r)); end
            weakLabels{end+1} = sprintf('%s%s (slack %.1e, lam %.1e)', ...
                rowLabel(r), nd, slackRow(r), R.lam_g(r)); %#ok<AGROW>
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
