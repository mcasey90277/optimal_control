function [rv0, rvf, meta] = insertion_states(target, criterion)
% INSERTION_STATES  Single source of truth for the GTO departure (rv0) and the
% tulip/ELFO insertion (rendezvous) state (rvf) used by every low-thrust
% pipeline. Declaring endpoints here (instead of threading them implicitly from
% a seed .mat) makes them explicit, changeable, and drift-checkable.
%
% INPUTS:
%   target    - 'tulip' | 'elfo' [char]
%   criterion - (optional) tulip: 'campaign'(default)|'maxydot'|'apoapsis'
%                          elfo:  'nearest'(default)|'apolune'|'perilune' [char]
%
% OUTPUTS:
%   rv0  - GTO departure state, ND rotating frame [1x6]
%   rvf  - insertion (rendezvous) state, ND rotating frame [1x6]
%   meta - struct: .target .criterion .label (label for filenames/provenance)
%
% NOTE: 'campaign'/'nearest' reproduce EXACTLY what the current seeds hold
% (Option A). Using an alternate criterion requires a matching energy seed --
% the consumer drivers' drift guard will fail loudly until one exists.
%
% REFERENCES:
%   [1] gto_tulip_endpoints.m (max-ydot tulip point + trace);
%   [2] gto_elfo_endpoints.m  (ELFO apolune/perilune/nearest).

if nargin < 2 || isempty(criterion)
    switch lower(target)
        case 'tulip', criterion = 'campaign';
        case 'elfo',  criterion = 'nearest';
        otherwise, error('insertion_states:target','unknown target %s', target);
    end
end
p = cr3bp_lt_params(25e-3, 15, 2100);
muStar = p.muStar;  tStar = p.tStar;  lStar = p.lStar; %#ok<NASGU>

% --- GTO departure (shared by all pipelines) --------------------------------
rv0 = [0.00349629072294633, -0.0072962582600817, 0, ...
       4.19147893803368, 8.98865558978329, 0];          % GTO departure (exact)
% -- to regenerate rv0 from the GTO orbital elements, uncomment: -------------
% muEarth = 6.67384e-20*(1-muStar)*(5.9736e24 + 7.35e22);
% sma = (6378+350 + 6378+35786)/2;  ecc = (35786-350)/(2*sma);
% [r0,v0] = pumpkyn.cr3bp.orb2eci(muEarth, [sma,ecc,0,-25*pi/180,0,0], 2);
% rv0 = pumpkyn.cr3bp.fromPCI(0, [r0,v0], muStar, tStar, lStar, 1);

here = fileparts(mfilename('fullpath'));

% --- insertion (rendezvous) state -------------------------------------------
switch lower(target)
  case 'tulip'
    switch lower(criterion)
      case 'campaign'
        rvf = [1.00658107295709, 0.0425745746906059, -0.0557780910480905, ...
               -0.16004281347248, 0.0665702939657711, -0.260455693516549];
        label = 'tulipCampaign';
      case 'maxydot'
        addpath(here);  [~, rvf] = gto_tulip_endpoints(p);      % max-ydot
        label = 'tulipMaxYdot';
      case 'apoapsis'
        addpath(here);  [~, ~, tr] = gto_tulip_endpoints(p);
        [~, idx] = min(vecnorm(tr(:,4:6), 2, 2));               % slowest point
        rvf = tr(idx, 1:6);
        label = 'tulipApoapsis';
      otherwise, error('insertion_states:crit','unknown tulip criterion %s', criterion);
    end
  case 'elfo'
    addpath(fullfile(here, '..', 'elfo'));                      % gto_elfo_endpoints
    switch lower(criterion)
      case 'nearest'
        [~, rvfTul] = insertion_states('tulip', 'campaign');    % ref = tulip campaign
        [~, rvf] = gto_elfo_endpoints(p, struct('point','nearest','ref',rvfTul));
        label = 'elfoNearest';
      case 'apolune'
        [~, rvf] = gto_elfo_endpoints(p, struct('point','apolune'));
        label = 'elfoApolune';
      case 'perilune'
        [~, rvf] = gto_elfo_endpoints(p, struct('point','perilune'));
        label = 'elfoPerilune';
      otherwise, error('insertion_states:crit','unknown elfo criterion %s', criterion);
    end
  otherwise, error('insertion_states:target','unknown target %s', target);
end
rvf  = rvf(:).';
meta = struct('target',lower(target),'criterion',lower(criterion),'label',label);
end
