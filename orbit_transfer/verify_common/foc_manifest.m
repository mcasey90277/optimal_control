function man = foc_manifest(name)
% FOC_MANIFEST  Campaign-level state and control dimension registry.
%
% Struct lookup for orbit-transfer campaign properties: state dimension (nx),
% control dimension (nu), row indices for direction/thrust/mass/time within
% the state vector, horizon type (fixed-tf or free-tf with scaling), and
% autonomy (whether the OCP has explicit time-dependence). Empty [] fields
% mean "not applicable" for that campaign. Used by the generic first-order
% optimality gate to route verification logic per campaign.
%
% INPUTS:
%   name - campaign identifier [char]; valid values: 'earth_mee', 'earth_cr3bp',
%          'tulip', 'elfo_fuel', 'elfo_mintime', 'toy'
%
% OUTPUTS:
%   man  - struct with fields [struct]:
%          .name              [char] campaign name
%          .nx                [int] state dimension
%          .nu                [int] control dimension
%          .dirRows           [1x3 or []] indices of X, Y, Z in state vector
%          .thrRow            [int or []] row index of thrust magnitude (or [] if control-indexed)
%          .massRow           [int or []] row index of mass in state
%          .timeRow           [int or []] row index of time in state
%          .autonomous        [logical] true if system is time-autonomous.
%                             DECLARATIVE ONLY (final-review fix I3, 2026-07-25):
%                             no code in foc_check.m/foc_report.m currently
%                             reads this field -- it is metadata recorded for
%                             a human/future consumer, not yet wired to any
%                             gate. The time-costate line foc_report prints is
%                             INFORMATIONAL for every campaign regardless of
%                             this flag (no PASS/FAIL is computed from it for
%                             any campaign, autonomous or not).
%          .horizonKind       [char] 'fixedtf' | 'freetf-cscale' | 'none'
%          .massFreeAtTf      [logical] true if final mass is a free variable
%
% REFERENCES:
%   [1] casadi_lt_mee.m header state/control layout (2026-07-25)
%   [2] casadi_minfuel_sundman.m header state/control layout (2026-07-25)
%   [3] casadi_energy_freetf.m header state/control layout (2026-07-25)
%   [4] casadi_mintime_freetf.m header state/control layout (2026-07-25)

base = struct('dirRows',[1 2 3], 'thrRow',4, 'autonomous',true, ...
              'horizonKind','fixedtf', 'massFreeAtTf',true);
switch lower(name)
    case 'earth_mee'
        man = base; man.name='earth_mee'; man.nx=7; man.nu=4; man.massRow=6; man.timeRow=7;
    case 'earth_cr3bp'
        man = foc_manifest('earth_mee'); man.name='earth_cr3bp'; man.autonomous=false;
    case 'tulip'
        man = base; man.name='tulip'; man.nx=8; man.nu=4; man.massRow=7; man.timeRow=8;
    case 'elfo_fuel'
        man = base; man.name='elfo_fuel'; man.nx=9; man.nu=4; man.massRow=7; man.timeRow=8;
        man.horizonKind = 'freetf-cscale';
    case 'elfo_mintime'
        man = foc_manifest('elfo_fuel'); man.name='elfo_mintime'; man.nu=3; man.thrRow=[];
    case 'toy'   % 1D double-integrator bang-bang used by test_foc_check_toy
        man = struct('name','toy','nx',2,'nu',1,'dirRows',[],'thrRow',1, ...
            'massRow',[],'timeRow',[],'autonomous',true,'horizonKind','none', ...
            'massFreeAtTf',false);
    otherwise, error('foc_manifest:unknown', 'unknown campaign %s', name);
end
end
