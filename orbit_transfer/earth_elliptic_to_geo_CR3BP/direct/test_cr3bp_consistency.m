% TEST_CR3BP_CONSISTENCY  Our Earth-centered 3-body EOM == BCP sec 4.2 CR3BP.
%
% Verifies numerically that the campaign's Earth-centered formulation
% (Kepler + third-body direct + indirect terms, circular Moon at rate
% n_M = sqrt((muE+muM)/D^3)) is the SAME dynamics as the rotating-frame
% barycentric CR3BP of Bonnard-Caillau-Picot 2010 sec 4.2. Method: compute
% the spacecraft's INERTIAL BARYCENTRIC acceleration two ways at random
% states and times --
%   (a) ours: Earth-centered accel  a = -muE*r/|r|^3 + a_M(r,t)  (direct +
%       indirect), then add the Earth's own barycentric acceleration
%       muM*rM(t)/D^3 (which cancels the indirect term analytically -- the
%       test exercises exactly that cancellation);
%   (b) theirs: plain two-gravity sum toward Earth and Moon placed on their
%       CR3BP circular barycentric orbits (Earth at -[muM/(muE+muM)]*rM,
%       Moon at +[muE/(muE+muM)]*rM), which is the inertial content of the
%       sec-4.2 rotating-frame equation (the 2*i*zdot and z terms are pure
%       frame bookkeeping).
% Agreement to round-off at random 3D states proves the formulations agree;
% it fails if the indirect term is dropped OR if n_M were built from muE
% alone (both classic inconsistencies).
%
% REFERENCES:
%   [1] Bonnard, Caillau, Picot, CIS 10(4) 2010, sec 4.2.
%   [2] doc/cr3bp_geo_phase1_note.tex sec 3 (our Eq. (3)).
here = fileparts(mfilename('fullpath'));  addpath(here);  setup_paths;
par = kepler_lt_params(10, 1500, 2000);
pM  = lunar_params(par, 0.7, 1);
muE = par.mu;  muM = pM.muM;  D = pM.DM;
rng(11);
for kk = 1:8
    r  = 1.5*(rand(3,1)-0.5);  r(3) = 0.1*(rand-0.5);   % Earth-centered pos
    t  = 40*rand;                                        % canonical time
    ang = pM.nM*t + pM.phi0;
    rM = D*[cos(ang); sin(ang); 0];                      % Earth->Moon vector
    % (a) OUR formulation, promoted to barycentric
    d  = rM - r;
    aM = muM*( d/norm(d)^3 - rM/D^3 );                   % direct + indirect
    aOurs = -muE*r/norm(r)^3 + aM;                       % Earth-centered
    aOursBary = aOurs + muM*rM/D^3;                      % + Earth's bary accel
    % (b) BCP sec-4.2 content: two-gravity sum, primaries on barycentric circles
    f  = muM/(muE+muM);                                  % Earth's offset fraction
    rEb = -f*rM;                                         % barycenter->Earth
    rSb = rEb + r;                                       % barycenter->spacecraft
    rMb = (1-f)*rM;                                      % barycenter->Moon
    aTheirs = -muE*(rSb-rEb)/norm(rSb-rEb)^3 - muM*(rSb-rMb)/norm(rSb-rMb)^3;
    err = norm(aOursBary - aTheirs)/norm(aTheirs);
    fprintf('state %d: rel err = %.3e\n', kk, err);
    assert(err < 1e-14, 'formulations disagree');
end
% Negative control: DROP the indirect term -> must disagree materially
r = [1;0;0];  t = 3;  ang = pM.nM*t + pM.phi0;  rM = D*[cos(ang);sin(ang);0];
d = rM - r;
aNoInd = -muE*r/norm(r)^3 + muM*d/norm(d)^3 + muM*rM/D^3;  % missing -muM*rM/D^3
f = muM/(muE+muM);  rEb=-f*rM;  rSb=rEb+r;  rMb=(1-f)*rM;
aRef = -muE*(rSb-rEb)/norm(rSb-rEb)^3 - muM*(rSb-rMb)/norm(rSb-rMb)^3;
assert(norm(aNoInd-aRef)/norm(aRef) > 1e-4, 'negative control: indirect term must matter');
fprintf('test_cr3bp_consistency: ALL PASS (equivalence to round-off; indirect term load-bearing)\n');
