function GH = highgamma_select(resDir)
% HIGHGAMMA_SELECT  Pick, for every high-gamma record of the three-arm race
% (FINDINGS 28), the arm to package: the SHARPEST solution among the arms
% that reached the floor. Sharpness = ramp width of the throttle law at the
% arm's deepest rung -- huber 0 (a jump), huberc delta_floor, eps 2p -- so
% the entry closest to the bang-bang limit wins; ties fall to eps (the
% catalog's original family). Records where no arm reached the floor are
% NOT selected (they are not min-fuel solutions).
%
% Returns records in the minfuel_grid.mat G-format plus .family/.delta so
% build_minfuel_catalog and run_conj_fixedtf_sweep treat them like grid
% records (source 'highgamma').
%
% INPUTS:
%   resDir - direct/results folder [char] (default: this campaign's)
%
% OUTPUTS:
%   GH - struct array: .cellIdx .gamma .pDeepest .delta .mfFuel .mfEnergy
%        .dMfPct .coastFrac .z .Y .family .width .source ('highgamma')
%        .altArms (cellstr of the other arms that also arrived)
%
% REFERENCES:
%   [1] run_highgamma_race.m (the arms), FINDINGS 28 (the table).

if nargin < 1
    resDir = fullfile(fileparts(mfilename('fullpath')), 'direct', 'results');
end
Lh = load(fullfile(resDir, 'highgamma_race.mat'));  H = Lh.H;
pFloor = 0.0015;
keys = unique(arrayfun(@(h) sprintf('%d_%d_%.10f', h.cellIdx, h.gamma), H, 'UniformOutput', false));
GH = struct([]);
for k = 1:numel(keys)
    rows = H(strcmp(arrayfun(@(h) sprintf('%d_%d_%.10f', h.cellIdx, h.gamma), H, 'UniformOutput', false), keys{k}));
    arrived = rows(arrayfun(@(h) isfinite(h.pDeepest) && h.pDeepest <= pFloor, rows));
    if isempty(arrived), continue, end
    w = arrayfun(@(h) rampWidth(h), arrived);
    [~, ib] = min(w + 1e-12*strcmp({arrived.family}, 'huber'));   % ties -> not huber (eps preferred)
    tie = find(abs(w - w(ib)) < 1e-15);
    if numel(tie) > 1 && any(strcmp({arrived(tie).family}, 'eps')), ib = tie(strcmp({arrived(tie).family}, 'eps')); ib = ib(1); end
    b = arrived(ib);
    tag = sprintf('c%d%d_g%03d_%s', b.cellIdx, round(100*b.gamma), b.family);
    A = load(fullfile(resDir, ['minfuel_hg_' tag '.mat']));  A = A.out.arms.(b.family);
    g = struct('cellIdx', b.cellIdx, 'gamma', b.gamma, 'pDeepest', A.p(end), 'delta', A.delta(end), ...
        'mfFuel', A.mf(end), 'mfEnergy', b.mfEnergy, 'dMfPct', 100*(A.mf(end) - b.mfEnergy), ...
        'coastFrac', A.coastFrac(end), 'z', A.z, 'Y', A.Y{end}, 'family', b.family, ...
        'width', w(ib), 'source', 'highgamma', 'altArms', {setdiff({arrived.family}, b.family)});
    if isempty(GH), GH = g; else, GH(end+1) = g; end %#ok<AGROW>
end
end

function w = rampWidth(h)
% RAMPWIDTH  Throttle-law ramp width at the arm's deepest rung.
% INPUTS: h race row.  OUTPUTS: w double.
switch h.family
    case 'huber',  w = 0;
    case 'huberc', w = h.delta;
    case 'eps',    w = 2*h.pDeepest;
    otherwise,     w = Inf;
end
end
