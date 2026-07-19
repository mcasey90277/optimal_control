function v = sosc_decide(K, AS, IN)
% SOSC_DECIDE  Map KKT/active-set/inertia results to a verdict (DESIGN sec 5).
% OUTPUTS: v - struct .verdict .reason .status
if ~K.pass || ~K.signOK
    v.verdict='ERROR'; v.reason='KKT residual/sign check failed';
elseif ~AS.licq || AS.nWeak>0 || IN.nzero>0
    parts={}; if ~AS.licq, parts{end+1}='LICQ fails'; end
    if AS.nWeak>0, parts{end+1}=sprintf('%d weakly-active: %s',AS.nWeak,strjoin(AS.weakLabels,'; ')); end
    if IN.nzero>0, parts{end+1}=sprintf('%d zero KKT eigenvalue(s)',IN.nzero); end
    v.verdict='INCONCLUSIVE'; v.reason=strjoin(parts,'; ');
elseif IN.subspaceOK
    v.verdict='PASS'; v.reason='reduced Hessian PD on the critical cone (strict local min)';
else
    v.verdict='FAIL'; v.reason=sprintf('indefinite reduced Hessian: inertia [%d %d %d] != expected [%d %d %d]', ...
        IN.npos,IN.nneg,IN.nzero,IN.expected(1),IN.expected(2),IN.expected(3));
end
switch v.verdict
    case 'PASS',         v.status='certified-sosc';
    case 'FAIL',         v.status='feasible-only';
    otherwise,           v.status='certified-feasibility+sosc-inconclusive';
end
end
