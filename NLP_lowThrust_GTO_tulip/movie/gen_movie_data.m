% gen_movie_data  Produce the trajectory data for the min-fuel solution movie.
%   - min-time full GTO->tulip spiral (faint context backdrop)
%   - the verified min-fuel NLP arrival-leg solution (bold, animated)
%   - a tulip reference arc
% Saves everything to minfuel_movie_data.mat.

addpath('/Users/msc/Desktop/optimal_control/lowThrust_GTO_tulip');
addpath('/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip');
run('/Users/msc/Desktop/optimal_control/NLP_lowThrust_GTO_tulip/setup_paths.m');  % pumpkyn on path

muStar = 0.012150585609624;
lStar  = 389703.264829278;
tStar  = 382981.289129055;
m0kg   = 15;
g0     = 9.80665*tStar^2/(1000*lStar);
Tmax   = (0.025/m0kg)*tStar^2/(lStar*1000);
c      = (2100/tStar)*g0;

muEarth = 6.67384e-20*(1 - muStar)*(5.9736E24 + 7.35E22);
sma = (6378+350 + 6378+35786)/2;
ecc = (35786 - 350)/(2*sma);
[r0d, v0d] = pumpkyn.cr3bp.orb2eci(muEarth, [sma, ecc, 0, -25*pi/180, 0, 0], 2);
rv0Full = pumpkyn.cr3bp.fromPCI(0, [r0d, v0d], muStar, tStar, lStar, 1);
[~, x0Tulip] = pumpkyn.cr3bp.getTulip((5/6)*2*pi, 7, -1, 1e-12);

tfMinFull = 6.2906939607;
zMinTime  = [190.4760481; -79.7060409; -0.4298691037; 0.3011592775; ...
              0.5866700046; -0.007117348902; 4.329378839];
optsInt = odeset('RelTol', 1e-12, 'AbsTol', 1e-14);

% --- min-time full spiral (backdrop) ---------------------------------------
fprintf('integrating min-time spiral...\n');
[tauF, yF] = ode113(@lt_pmp_eom, linspace(0, tfMinFull, 5000), ...
                    [rv0Full(:); 1; zMinTime], optsInt, Tmax, c, muStar);
legStart = 4.0;
[~, iLeg] = min(abs(tauF - legStart));   % index where the min-fuel leg begins

% --- tulip reference arc ---------------------------------------------------
tYul = []; yTul = [];
try
    [tYul, yTul] = pumpkyn.cr3bp.prop(2*pi, x0Tulip, muStar);
catch mePropErr
    fprintf('tulip reference prop failed: %s\n', mePropErr.message);
end

% --- verified min-fuel NLP arrival leg -------------------------------------
fprintf('solving min-fuel NLP leg (N=3000)...\n');
out = NLP_lowThrust_GTO_Tulip_minfuel(3000, 1.3, 4.0, false);

outFile = fullfile(fileparts(mfilename('fullpath')), 'minfuel_movie_data.mat');
save(outFile, 'tauF', 'yF', 'iLeg', 'legStart', 'out', 'x0Tulip', ...
     'tYul', 'yTul', 'muStar', 'lStar', 'tStar', 'Tmax', 'c', 'm0kg', 'tfMinFull');
fprintf('SAVED %s\n', outFile);
fprintf('leg: mProp=%.4f kg, flag=%d, burnFrac=%.1f%%, maxDefect=%.2g\n', ...
        out.mProp_kg, out.exitflag, 100*out.burnFrac, out.maxDefect);
