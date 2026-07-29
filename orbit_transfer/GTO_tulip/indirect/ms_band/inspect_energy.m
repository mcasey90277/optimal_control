% INSPECT_ENERGY  Layout check of the Branch-B candidate seed file.
S = load('/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/sundman_minfuel/results/energy/energy_f1120.mat');
fn = fieldnames(S);
for k = 1:numel(fn)
    v = S.(fn{k});
    fprintf('%-12s  %s  %s\n', fn{k}, mat2str(size(v)), class(v));
end
if isfield(S, 'out')
    fprintf('--- out fields ---\n');
    fo = fieldnames(S.out);
    for k = 1:numel(fo)
        v = S.out.(fo{k});
        fprintf('  %-14s  %s  %s\n', fo{k}, mat2str(size(v)), class(v));
    end
    if isfield(S.out, 'X')
        fprintf('X(8,end)=%.6f  factor=%.4f\n', S.out.X(8,end), S.factor);
        fprintf('U4 range [%.3g %.3g] mean %.3g\n', min(S.out.U(4,:)), max(S.out.U(4,:)), mean(S.out.U(4,:)));
    end
    if isfield(S.out, 'lamDef'), fprintf('lamDef size %s\n', mat2str(size(S.out.lamDef))); end
    if isfield(S.out, 'dV'), fprintf('dV=%.4f switches=%s\n', S.out.dV, mat2str(S.out.switches)); end
    if isfield(S.out, 'eps'), fprintf('out.eps=%.4g\n', S.out.eps); end
end
if isfield(S, 'eps'), fprintf('top eps=%.4g\n', S.eps); end
