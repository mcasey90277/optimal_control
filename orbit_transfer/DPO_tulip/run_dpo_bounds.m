function B = run_dpo_bounds()
%% Purpose:
%
%   Admissibility survey for the DPO -> TULIP campaign: finds the
%   "reasonable" distant prograde orbits by Darin's criteria (periselene
%   altitude >= 500 km, whole orbit within 100 Mm of the Moon), using the
%   shared costate_common/survey_family_bounds engine.
%
%   Candidates are pumpkynPie's NATIVE catalog members -- getDPO([]) returns
%   every seed period in its lookup table. The halo survey taught us that
%   interpolated members can fail cont_np and produce garbage metrics; the
%   survey's periodicity guard excludes them, and native seeds avoid most of
%   the problem outright.
%
%% Inputs:
%
%   none (run from MATLAB after pumpkynPie startup, or via the shell line
%   in the header of run_dpo_catalog.m)
%
%% Outputs:
%
%  B                        struct                  survey_family_bounds
%                                                   output, also saved to
%                                                   direct/results/
%                                                   dpo_bounds.mat
%
%% Revision History:
%  M. Casey                                                   (c) 08/07/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

here = fileparts(mfilename('fullpath'));
addpath(fullfile(fileparts(here), 'costate_common'));
resDir = fullfile(here, 'direct', 'results');
if ~isfolder(resDir), mkdir(resDir); end

%% Native catalog members (tau vector straight from the lookup table):
tauAll = pumpkynPie.cr3bp.getDPO([]);
tauAll = tauAll(:).';
fprintf('getDPO catalog: %d native members, tau %.3f-%.3f ND\n', ...
        numel(tauAll), min(tauAll), max(tauAll));

%% Survey every native member (props are cheap; no need to subsample):
grid_ = struct('tau', num2cell(tauAll));
B = survey_family_bounds('dpo', grid_, fullfile(resDir, 'dpo_bounds.mat'));

%% Admissible box summary, in the units the front door wants:
if any(B.admissible)
    tAdm = arrayfun(@(r) r.params.tau, B.rows(B.admissible));
    fprintf(['\nADMISSIBLE DPO BOX: tau %.3f-%.3f ND ', ...
             '(%.2f-%.2f days), %d members\n'], ...
            min(tAdm), max(tAdm), min(tAdm)*4.432654, ...
            max(tAdm)*4.432654, numel(tAdm));
end
end
