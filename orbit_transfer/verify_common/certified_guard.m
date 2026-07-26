function info = certified_guard(res, spec)
% CERTIFIED_GUARD  Refuse to verify against a re-solve that is not the certified point.
%
% Every FOC/PMP driver in this repo verifies a SAVED solution by warm re-solving
% it to recover multipliers. That re-solve can land somewhere other than the
% certified point, and gating optimality conditions on the wrong point produces a
% confident, wrong certificate. This is the single gate that stands between those
% two outcomes, in three parts:
%
%   1. SUCCESS CLASS - IPOPT actually converged. An unconverged point has no
%      usable multipliers at all.
%   2. FEASIBILITY   - maxDefect is machine-tight. Multipliers at an infeasible
%      point are not a KKT certificate.
%   3. CERTIFIED QUANTITY, ONE-SIDED - the objective that defines the row did not
%      get WORSE. One-sided is the substance: a warm re-solve that finds a nearby
%      BETTER optimum is legitimate to verify (it is a machine-tight extremal of
%      the SAME problem -- same thrust, same pinned horizon, same endpoints) and
%      is flagged via info.improved rather than refused. Only degradation is
%      refused. Both directions are supported because the campaigns differ:
%      min-fuel maximizes m_f (higher is better), min-time minimizes t_f (lower
%      is better).
%
% NOT a node-wise drift check. Two machine-tight extrema of the same problem can
% differ node-by-node while both being valid; requiring node agreement rejects
% almost every honest re-solve. The certified QUANTITY is the invariant.
%
% WHY THIS IS SHARED. It was duplicated across four drivers in four campaigns
% (run_foc_tulip, run_foc_elfo x2, refresh_duals_mee, refresh_duals_cr3bp) with
% identical thresholds and identical logic. The apparent divergence between them
% -- two testing out.success alone, two also whitelisting the IPOPT status -- was
% NOT a real difference: casadi_lt_mee folds the status whitelist into
% out.success itself, while the Sundman solvers set success=true on any
% non-throwing solve, so their callers had to apply it. This function applies the
% whitelist unconditionally, which is a no-op for the former and required for the
% latter.
%
% INPUTS:
%   res  - re-solve result [struct]:
%          .success     IPOPT reported success                     [logical]
%          .ipoptStatus solver return status                       [char]
%          .maxDefect   max collocation defect                     [scalar]
%          .value       the RESOLVED certified quantity            [scalar]
%   spec - guard specification [struct]:
%          .caller    driver name; prefixes every error id         [char]
%          .label     row/file tag used in messages                [char]
%          .saved     the CERTIFIED value of the quantity          [scalar]
%          .name      quantity name for messages, e.g. 'm_f'       [char]
%          .errName   quantity token in the error id, e.g. 'mass'  [char]
%          .better    'higher' (min-fuel m_f) | 'lower' (min-time t_f) [char]
%          .feasTol   defect tolerance                    [scalar, default 1e-8]
%          .tol       relative tolerance on the quantity  [scalar, default 1e-6]
%          .units     optional unit suffix for messages   [char, default '']
%
% OUTPUTS:
%   info - [struct] .ipoptStatus .maxDefect .saved .resolved .relChange
%          .improved (true if the re-solve found a strictly BETTER point)
%
% THROWS:
%   <caller>:resolveFailed    not a success-class IPOPT status
%   <caller>:notFeasible      maxDefect above feasTol
%   <caller>:<errName>Degraded  certified quantity moved the wrong way
%
% REFERENCES:
%   [1] earth_elliptic_to_geo/direct/reproduce/verify_row.m -- origin of the
%       one-sided convention (minimum fuel means MAXIMIZE final mass, so a
%       higher mass can never be a failure).
%   [2] verify_common/tests/test_certified_guard.m -- pins all three refusal
%       branches, both directions, and the improved/unchanged paths.

if ~isfield(spec, 'feasTol') || isempty(spec.feasTol), spec.feasTol = 1e-8;  end
if ~isfield(spec, 'tol')     || isempty(spec.tol),     spec.tol     = 1e-6;  end
if ~isfield(spec, 'units')   || isempty(spec.units),   spec.units   = '';    end

% A mislabelled direction would silently invert the whole gate -- refusing good
% points and accepting degraded ones -- so it is validated, not defaulted.
if ~any(strcmp(spec.better, {'higher', 'lower'}))
    error('certified_guard:badSpec', ...
        'spec.better must be ''higher'' or ''lower'', got ''%s''', spec.better);
end

lbl = spec.label;
if ~isempty(lbl), lbl = [lbl ': ']; end

% --- 1. success class -----------------------------------------------------
successClass = res.success && any(strcmp(res.ipoptStatus, ...
    {'Solve_Succeeded', 'Solved_To_Acceptable_Level'}));
if ~successClass
    error([spec.caller ':resolveFailed'], ...
        ['%swarm re-solve did not reach a success-class IPOPT status (%s) -- ' ...
         'cannot gate the FOC/KKT check on an unconverged point. The certified ' ...
         'primal is unaffected; only verification is blocked for this row.'], ...
        lbl, res.ipoptStatus);
end

% --- 2. feasibility -------------------------------------------------------
if res.maxDefect > spec.feasTol
    error([spec.caller ':notFeasible'], ...
        ['%swarm re-solve maxDefect %.3e > feasTol %.3e -- the recovered point ' ...
         'is not machine-tight, so its multipliers are not a usable KKT ' ...
         'certificate.'], lbl, res.maxDefect, spec.feasTol);
end

% --- 3. certified quantity, one-sided ------------------------------------
saved = spec.saved;  got = res.value;
scale = max(abs(saved), 1);              % absolute floor: saved may be ~0
relChange = abs(got - saved) / max(abs(saved), 1e-30);

info = struct('ipoptStatus', res.ipoptStatus, 'maxDefect', res.maxDefect, ...
              'saved', saved, 'resolved', got, 'relChange', relChange, ...
              'improved', false);

if strcmp(spec.better, 'higher')
    degraded = got < saved - spec.tol*scale;
    improved = got - saved >  spec.tol*scale;
else
    degraded = got > saved + spec.tol*scale;
    improved = saved - got >  spec.tol*scale;
end

if degraded
    error([spec.caller ':' spec.errName 'Degraded'], ...
        ['%swarm re-solve %s DEGRADED (%.6f -> %.6f%s, %.3e relative) -- the ' ...
         'recovered point is worse than the certified one; refusing to verify ' ...
         'against it.'], lbl, spec.name, saved, got, spec.units, relChange);
elseif improved
    info.improved = true;
    fprintf(['  [%s] NOTE %swarm re-solve improved %s (%.6f -> %.6f%s, %.3e ' ...
             'relative). Verifying the IMPROVED extremal; the gates below ' ...
             'describe that point, not the published row.\n'], ...
            spec.caller, lbl, spec.name, saved, got, spec.units, relChange);
end
end
