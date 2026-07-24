% PHASE_CONTROL_SENSITIVITY  Quantify control-law differences across the
% phi0 sweep (10 N): switch-time shifts, throttle deltas, structure check.
%
% Loads the four phi0-swept certified products and compares each control law
% against the phi0=0 reference: bang-bang switch count (structure), switch
% CROSSING TIMES (linear interpolation of the throttle through 0.5 on the
% time grid, reported in minutes), and nodal throttle differences. Feeds the
% "control-law sensitivity" paragraph of the note's lunar-phase section:
% same 19-switch structure everywhere, switch shifts <= ~1.3 min on a 5.3 d
% transfer, zero burn/coast node flips, and the phi0=pi law nearly
% coinciding with phi0=0 (the quadrupole symmetry at the control level).
%
% INPUTS:  none        OUTPUTS: none (table printed)
% REFERENCES:
%   [1] doc/cr3bp_geo_phase1_note.tex, sec 'Lunar-phase dependence'
%       (control-law sensitivity paragraph -- the numbers printed here).
here   = fileparts(mfilename('fullpath'));
resDir = fullfile(here, 'results');
runs = {'phi0=0     ','cr3bp_T10N_phi0_fuel.mat'; 'phi0=pi/2  ','cr3bp_T10N_phiPi2_fuel.mat'; ...
        'phi0=pi    ','cr3bp_T10N_phiPi_fuel.mat'; 'phi0=3pi/2 ','cr3bp_T10N_phi3Pi2_fuel.mat'};
S0 = load(fullfile(resDir, runs{1,2}));
t0 = local_cross(S0.t_days, S0.throttle);
fprintf('reference %s: %d switches\n', strtrim(runs{1,1}), numel(t0));
for k = 2:4
    S  = load(fullfile(resDir, runs{k,2}));
    tk = local_cross(S.t_days, S.throttle);
    fprintf('%s: %d switches; ', runs{k,1}, numel(tk));
    if numel(tk) == numel(t0)
        d = (tk - t0) * 24 * 60;
        fprintf('shift vs ref: max %+.1f min, mean|.| %.1f min; ', ...
            d(find(abs(d)==max(abs(d)),1)), mean(abs(d)));
    else
        fprintf('STRUCTURE CHANGED; ');
    end
    du = abs(S.throttle - S0.throttle);
    fprintf('max|dthr|=%.3f, flips=%.2f%%\n', max(du), 100*mean(du>0.5));
end

function tc = local_cross(t, u)
% Throttle 0.5-crossing times by linear interpolation [days].
s = u > 0.5;  ix = find(diff(s) ~= 0);
tc = zeros(numel(ix),1);
for m = 1:numel(ix)
    a = ix(m);  tc(m) = t(a) + (0.5-u(a))*(t(a+1)-t(a))/(u(a+1)-u(a));
end
end
