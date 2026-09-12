%% COMPARE_RESWEEP  Rebuilt (results_resweep) vs shipped (results) 70 mN catalog.
%   Entries and costates must match bitwise (the sheet stage is deterministic
%   and the ribs are reused); what is EXPECTED to change is the second-order
%   sheet: the 42 old near-miss cells re-measured by the resolved scan.
base = '/Users/msc/Desktop/optimal_control/orbit_transfer/DRO_tulip/indirect/results';
rr   = [base '_resweep'];
A = load(fullfile(base, 'costate_catalog_dro_tulip_70mN.mat'));  fa = fieldnames(A);  ca = A.(fa{1});
B = load(fullfile(rr,   'costate_catalog_dro_tulip_70mN.mat'));  fb = fieldnames(B);  cb = B.(fb{1});
sa = ca.sheets(1);  sb = cb.sheets(1);
fprintf('entries: shipped %d, rebuilt %d\n', ca.n_entries, cb.n_entries);
onlyA = find(sa.has_solution & ~sb.has_solution);  onlyB = find(~sa.has_solution & sb.has_solution);
fprintf('cells only in shipped: %d, only in rebuilt: %d\n', numel(onlyA), numel(onlyB));
both = find(sa.has_solution & sb.has_solution);
dz = zeros(1, numel(both));
for k = 1:numel(both)
    dz(k) = max(abs(sa.z8(:, sa.entry_index(both(k))) - sb.z8(:, sb.entry_index(both(k)))));
end
fprintf('common cells %d: max |dz8| %.2e\n', numel(both), max(dz));
fprintf('conj_pass shipped %d / rebuilt %d on common cells\n', nnz(sa.conj_pass(both) == 1), nnz(sb.conj_pass(both) == 1));
Ua = load(fullfile(base, 'audit_70mN.mat'));  Ub = load(fullfile(rr, 'audit_70mN.mat'));
fprintf('audit: shipped %d ok / %d bad; rebuilt %d ok / %d bad\n', Ua.A.nOk, Ua.A.nBad, Ub.A.nOk, Ub.A.nBad);
if Ub.A.nBad > 0, fprintf('   %s\n', Ub.A.problems{:}); end
Sa = load(fullfile(base, 'arrival_sheet_70mN.mat'));  Sb = load(fullfile(rr, 'arrival_sheet_70mN.mat'));
fprintf('sheet spine: certified columns shipped %d / rebuilt %d\n', nnz(isfinite(Sa.S.TF)), nnz(isfinite(Sb.S.TF)));
fprintf('certified per column shipped [%s] rebuilt [%s]\n', num2str(Sa.S.nCert(:)'), num2str(Sb.S.nCert(:)'));
% refusals under the widened stack, by reason
for j = 1:numel(Sb.S.cand)
    for k = 1:numel(Sb.S.cand{j})
        c = Sb.S.cand{j}(k);
        if ~c.ok && (numel(Sa.S.cand{j}) >= k) && Sa.S.cand{j}(k).ok
            fprintf('   NEWLY REFUSED col %d cand %d (%.2f d): %s\n', j, k, c.tfDays, c.reason);
        end
    end
end
% the second-order sheet, before and after
f2 = {'conj_interior', 'conj_zero', 'conj_near_miss', 'conj_interior_cand', 'conj_multiplicity'};
for k = 1:numel(f2)
    va = sa.(f2{k})(both);  vb = sb.(f2{k})(both);
    fprintf('%-20s cells > 0: shipped %d, rebuilt %d\n', f2{k}, nnz(va > 0), nnz(vb > 0));
end
if isfield(sb, 'conj_unresolved')
    vu = sb.conj_unresolved(both);
    fprintf('%-20s cells > 0: rebuilt %d (NaN %d)\n', 'conj_unresolved', nnz(vu > 0), nnz(isnan(vu)));
    if any(vu > 0)
        [i1, i2] = ind2sub(size(sb.has_solution(:,:,1)), both(vu > 0));
        for k = 1:numel(i1), fprintf('   UNRESOLVED at cell (%d,%d)\n', i1(k), i2(k)); end
    end
end
fprintf('H6 worst: shipped %.2fx, rebuilt %.2fx; lift worst: shipped %.0fx, rebuilt %.0fx\n', ...
        min(sa.h6_margin(both)), min(sb.h6_margin(both)), min(sa.lift_margin(both)), min(sb.lift_margin(both)));
