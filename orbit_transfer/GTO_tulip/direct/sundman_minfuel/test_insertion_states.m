% TEST_INSERTION_STATES  Verify the insertion-point helper: the default criteria
% reproduce exactly what the existing seeds hold (so the drift guards pass with
% zero re-solve), and the alternate criteria return valid 6-states.
here = fileparts(mfilename('fullpath'));  cd(here);  setup_paths();  addpath('../elfo');
E  = load('results/energy/energy_f1120.mat');           % a tulip backbone
Ee = load('../elfo/results/energy_elfo_freetf.mat');    % the ELFO seed

[rv0,rvfC,mC] = insertion_states('tulip','campaign');
assert(norm(rv0  - E.rv0(:).') < 1e-12, 'rv0 != backbone rv0');
assert(norm(rvfC - E.rvf(:).') < 1e-12, 'tulip campaign rvf != backbone rvf');
assert(strcmp(mC.label,'tulipCampaign'), 'wrong tulip label');

[~,rvfN,mN] = insertion_states('elfo','nearest');
assert(norm(rvfN - Ee.rvf(:).') < 1e-12, 'elfo nearest rvf != ELFO seed rvf');
assert(strcmp(mN.label,'elfoNearest'), 'wrong elfo label');

[~,rvfM] = insertion_states('tulip','maxydot');
[~,rvfA] = insertion_states('tulip','apoapsis');
[~,rvfP] = insertion_states('elfo','apolune');
assert(all(isfinite(rvfM)) && numel(rvfM)==6, 'maxydot invalid');
assert(all(isfinite(rvfA)) && numel(rvfA)==6, 'apoapsis invalid');
assert(all(isfinite(rvfP)) && numel(rvfP)==6, 'apolune invalid');

% default criterion (omitted arg) == campaign / nearest
[~,rvfCd] = insertion_states('tulip');  assert(norm(rvfCd-rvfC)<1e-15,'tulip default');
[~,rvfNd] = insertion_states('elfo');   assert(norm(rvfNd-rvfN)<1e-15,'elfo default');
fprintf('TEST_INSERTION_STATES: PASS (defaults match seeds <1e-12; alternates valid)\n');
