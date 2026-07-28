function S = mesh_collect(L, opts)
% MESH_COLLECT  Assemble per-quantity series and orders from a mesh ladder.
%
% Turns a ladder (mesh_ladder_mee output) into the series the study reports:
% each tracked quantity at each refinement level, with an observed order where
% one is admissible.
%
% TWO CLASSIFICATIONS DECIDE WHAT MAY BE SAID ABOUT EACH QUANTITY, and both are
% carried in the output rather than left to the reader:
%
%   .class   'physical'   a property of the trajectory. May earn a convergence
%                         verdict.
%            'integer'    a count. Must STABILIZE; no order exists for it.
%            'diagnostic' a solver/discretization readout (KKT residual, dual
%                         statistics, sign percentages). Tracked and reported,
%                         but NEVER given a CONVERGED verdict alongside
%                         physical quantities -- these measure the NLP and its
%                         solution process, not the continuous OCP.
%
%   .richardsonOK  whether Richardson extrapolation is PERMITTED here. The
%                  study's policy forbids it for switch counts, raw switch
%                  times on an unaligned grid, raw endpoint defect multipliers,
%                  sign percentages, KKT residuals, and anything measured
%                  across a topology change. Richardson assumes a smooth
%                  asymptotic expansion in h; with switch locations moving
%                  between meshes the leading coefficient oscillates with mesh
%                  phase and a three-level extrapolation can be badly
%                  misleading. Where it is forbidden, .order.p and .order.rich
%                  are forced to NaN -- the series and its level-to-level
%                  deltas are still reported, because those are honest.
%
% SWITCH MATCHING USES A PER-SWITCH WINDOW, NOT A GLOBAL ONE. Each level's
% switches are matched against the FINEST level's, with switch i of level k
% given a window of 2x the physical length of ITS OWN bracketing interval on
% level k. This matters because these meshes are uniform in longitude or in a
% Sundman variable, never in time: the physical step was measured to vary 36-39x
% across the earth 10 N transfer. A global window built from the median step
% was measured to spuriously refuse a genuine pair whose switch had moved
% 0.382 while the median step was 0.156. Using level k's own local step (rather
% than always the coarsest level's) is the coarser of each compared pair, which
% is where the location uncertainty actually lives.
%
% INPUTS:
%   L    - ladder struct array from mesh_ladder_mee [1xL]. Levels with
%          .ok == false are DROPPED, with a warning: a failed level is not a
%          solution and must not enter a series.
%   opts - struct (optional):
%          .tolMult   per-switch window, in local steps   [default 2]
%          .m0kg      initial mass, for reporting m_f in kg when out.m_f_kg
%                     is absent [default: taken from out.m_f_kg]
%
% OUTPUTS:
%   S - struct:
%       .factors  [1xL] refinement multipliers actually used
%       .N        [1xL] interval counts
%       .nLevels  scalar
%       .hPhys    [1xL] MEDIAN physical step per level -- the study reports
%                 physical step sizes, not node-count multipliers, because
%                 "h = 1/N" is not the relevant h on a nonuniform mesh
%       .hPhysMax [1xL] max physical step per level
%       .q        struct of quantities; each .q.<name> is
%                 struct('vals',[1xL], 'order',<mesh_order out>, 'class',
%                 'richardsonOK', 'label', 'units')
%       .switchStability - struct:
%                 .counts        [1xL] switch count per level
%                 .stable        logical, all counts equal
%                 .matched       [1xL] pairs against the finest level
%                 .unmatchedSelf [1xL] switches on level k with no fine partner
%                 .unmatchedFine [1xL] fine switches with no level-k partner
%                 .maxAbsDt      [1xL] largest matched movement
%                 .medAbsDt      [1xL] median matched movement
%                 .tolMedian     [1xL] median window used
%                 .extraFine     {1xL} times of the unmatched fine switches
%       .switchTimes  {1xL} the per-level switch time vectors
%       .switchMeta   {1xL} the mesh_switch_times structs (method, fracs)
%
% REFERENCES:
%   [1] docs/superpowers/plans/2026-07-25-mesh-convergence-study.md, Task 4
%       and "Richardson policy (restricted)".
%   [2] mesh_switch_times.m, mesh_match_switches.m, mesh_order.m.

if nargin < 2, opts = struct(); end
tolMult = local_default(opts, 'tolMult', 2);

bad = ~[L.ok];
if any(bad)
    warning('mesh_collect:droppedLevels', ...
        ['dropping %d failed level(s) from the series -- a level that did not ' ...
         'converge is not a solution and must not be fitted'], sum(bad));
    L = L(~bad);
end
nL = numel(L);
assert(nL >= 2, 'mesh_collect:levels', ...
    'need at least 2 converged levels, got %d', nL);

S = struct();
S.factors = [L.factor];
S.N       = [L.N];
S.nLevels = nL;

% Physical step per level. Reported because a node-count multiplier is not h on
% a mesh that is uniform in longitude or in a Sundman variable.
S.hPhys    = nan(1,nL);
S.hPhysMax = nan(1,nL);
for k = 1:nL
    dt = diff(L(k).out.X(7,:));
    S.hPhys(k)    = median(dt);
    S.hPhysMax(k) = max(dt);
end

% --- the tracked quantities -------------------------------------------------
% Column meanings: name, extractor, class, richardsonOK, label, units.
spec = {
 'mf',               @(o,r) o.m_f_kg,          'physical',   true,  'final mass',              'kg'
 'dV_kms',           @(o,r) o.dV_kms,          'physical',   true,  'delta-V',                 'km/s'
 ... % termErr is a DIAGNOSTIC, not a physical observable, and the plan's
 ... % expectations table drops it explicitly: the terminal elements are
 ... % EQUALITY-constrained in the NLP, so this reads solver tolerance by
 ... % construction and can never measure discretization error. Measured on the
 ... % earth 10 N ladder it is 1.6e-35 -> 2.6e-36, which would have earned a
 ... % meaningless CONVERGED verdict had it stayed classed physical. The plan's
 ... % intended replacement -- terminal error of the INDEPENDENTLY PROPAGATED
 ... % control -- is a different quantity and is not implemented here.
 'termErr',          @(o,r) o.termErr,         'diagnostic', false, 'terminal elem err (=tol)', '-'
 'dL',               @(o,r) o.dL,              'physical',   true,  'longitude span',          'rad'
 'nSwitches',        @(o,r) o.switches,        'integer',    false, 'switch count',            '-'
 'maxDefect',        @(o,r) o.maxDefect,       'diagnostic', false, 'max defect',              '-'
 'kktStatInf',       @(o,r) r.kktStatInf,      'diagnostic', false, 'KKT stationarity',        '-'
 'lamMassEndMapped', @(o,r) r.lamMassEndMapped,'diagnostic', false, 'mapped transversality',   '-'
 'lamTimeCoV',       @(o,r) r.lamTimeCoV,      'diagnostic', false, 'time-costate CoV',        '-'
 'sdotMinRel',       @(o,r) r.sdotMinRel,      'diagnostic', false, 'min regular Sdot',        '-'
 'signPct',          @(o,r) r.signPct,         'diagnostic', false, 'dual sign law',           '%'
 'dirTanMax',        @(o,r) r.dirTanMax,       'diagnostic', false, 'primer misalign (max)',   'deg'
};

S.q = struct();
for row = 1:size(spec,1)
    name = spec{row,1};  get = spec{row,2};
    vals = nan(1,nL);
    for k = 1:nL
        try
            v = get(L(k).out, L(k).rep);
            if isempty(v), v = NaN; end
            vals(k) = double(v(1));
        catch
            vals(k) = NaN;                  % quantity absent on this level
        end
    end
    o = mesh_order(vals, S.factors);
    if ~spec{row,4}                          % Richardson forbidden here
        o.p = NaN;  o.rich = NaN;  o.pWindows = nan(size(o.pWindows));
        o.windowsConsistent = false;
    end
    S.q.(name) = struct('vals', vals, 'order', o, 'class', spec{row,3}, ...
        'richardsonOK', spec{row,4}, 'label', spec{row,5}, 'units', spec{row,6});
end

% --- switch structure -------------------------------------------------------
S.switchTimes = cell(1,nL);
S.switchMeta  = cell(1,nL);
locH          = cell(1,nL);
for k = 1:nL
    sw = local_switches(L(k));
    S.switchTimes{k} = sw.tSw;
    S.switchMeta{k}  = sw;
    t = L(k).out.X(7,:);
    if sw.n > 0
        locH{k} = t(sw.bracket(:,2)).' - t(sw.bracket(:,1)).';   % [1xk]
    else
        locH{k} = zeros(1,0);
    end
end

st = struct('counts', [], 'stable', false, 'matched', nan(1,nL), ...
    'unmatchedSelf', nan(1,nL), 'unmatchedFine', nan(1,nL), ...
    'maxAbsDt', nan(1,nL), 'medAbsDt', nan(1,nL), 'tolMedian', nan(1,nL), ...
    'extraFine', {cell(1,nL)});
st.counts = cellfun(@numel, S.switchTimes);
st.stable = all(st.counts == st.counts(1));

tFine = S.switchTimes{end};
for k = 1:nL
    if k == nL
        st.matched(k) = numel(tFine);  st.unmatchedSelf(k) = 0;
        st.unmatchedFine(k) = 0;  st.maxAbsDt(k) = 0;  st.medAbsDt(k) = 0;
        st.tolMedian(k) = NaN;  st.extraFine{k} = [];
        continue
    end
    if isempty(S.switchTimes{k}) || isempty(tFine)
        st.extraFine{k} = tFine;  continue
    end
    tolK = tolMult * max(locH{k}, realmin);        % per-switch window
    m = mesh_match_switches(S.switchTimes{k}, tFine, tolK);
    st.matched(k)       = m.n;
    st.unmatchedSelf(k) = numel(m.unmatchedA);
    st.unmatchedFine(k) = numel(m.unmatchedB);
    st.maxAbsDt(k)      = m.maxAbsDt;
    st.medAbsDt(k)      = median(abs(m.dt));
    st.tolMedian(k)     = median(tolK);
    st.extraFine{k}     = tFine(m.unmatchedB);
end
S.switchStability = st;

% Matched switch-time movement as its own series. Richardson is FORBIDDEN on
% it (raw switch times on an unaligned grid), but the series itself is the
% study's switch-time evidence, so it is carried with the others.
mv = st.maxAbsDt(1:end-1);
o  = struct('p',NaN,'rich',NaN,'dLast',NaN,'rel',NaN,'monotone',false, ...
            'pWindows',[],'windowsConsistent',false,'r',NaN);
if numel(mv) >= 2
    oo = mesh_order(mv, S.factors(1:numel(mv)));
    o.dLast = oo.dLast;  o.rel = oo.rel;  o.monotone = oo.monotone;  o.r = oo.r;
end
S.q.switchDrift = struct('vals', [mv NaN], 'order', o, 'class', 'physical', ...
    'richardsonOK', false, 'label', 'switch movement', 'units', 't');
end

% ---------------------------------------------------------------------------
function sw = local_switches(Lk)
% LOCAL_SWITCHES  Sub-grid switch times for one level, or an empty result.
% Uses the de-weighted switching function when the report carries one, and
% falls back to the throttle otherwise -- mesh_switch_times records which.
% INPUTS:  Lk - one ladder level    OUTPUTS: sw - mesh_switch_times struct
sg   = Lk.sigma(:).';
t    = Lk.out.X(7,:);
thr  = Lk.out.U(4,:);
burn = thr > 0.5;
Shat = [];
if isstruct(Lk.rep) && isfield(Lk.rep,'Sdeweighted'), Shat = Lk.rep.Sdeweighted; end
sw = mesh_switch_times(sg, t, burn, Shat, thr);
end

% ---------------------------------------------------------------------------
function v = local_default(s, f, dflt)
% LOCAL_DEFAULT  s.(f) if present and nonempty, else dflt.
if isstruct(s) && isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = dflt; end
end
