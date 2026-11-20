function ok = test_transfer_difficulty()
%% Purpose:
%
%   Unit test for transfer_difficulty, the scorecard of a thrust-ladder
%   sheet. Every expectation is hand-computed:
%
%     1. 4 x 4 x 3 sheet (rungs 3, 2, 1 N; OK pattern chosen so counts and
%        shares are known): nAccepted [15 14 8], shares, TF median / min /
%        max / spread of rung 1 (1.2, 1.1, 1.4, 0.25), wall median, the
%        ladder summary (lowestRungAllCells 3 at shareAll 0.9 and 2 at
%        0.85, lowestRungAnyCell 1, nFullLadders 8, nEmptyCells 1),
%        attMedian;
%     2. delta-V of a rung against an independent rocket-equation
%        bookkeeping in SI units (Isp g0 ln(m0/mf));
%     3. roughness: a constant sheet is 0; a 3 x 3 sheet with one spike
%        (TF 2 among 1) gives the median 0.25 (spike 0.5, four adjacent
%        cells 0.25, four others 0); a one-rung sheet (the engine's
%        dropped trailing singleton) is handled;
%     4. revolutions: a circular arc of 1.5 turns about the Moon in the
%        INERTIAL frame, written in the rotating frame, gives revMedian
%        1.5 to 1e-3, and about the Earth it does not (the primary matters);
%        a cell with DDONE false is excluded from nDirect and revMedian;
%        no sidecar gives NaN;
%     5. a sheet with no accepted cell returns NaN statistics without error;
%     6. conjPass: shareRefuted and nUndecided against counted values;
%     7. print: D.lines is a cell of char and opts.print = false is silent.
%
%% Outputs:
%
%   ok  logical  true when every check passes; also an error when false.
%
%% Revision History:
%
%   M. Casey (c) 11/19/2026
%
%   Copyright Coorbital Inc.
%

here = fileparts(mfilename('fullpath'));
addpath(fileparts(here));
tmp = tempname;  mkdir(tmp);
cleanup = onCleanup(@() rmdir(tmp, 's'));
ok = true;
op = struct('print', false);
muStar = 0.012150585609624;

% -------- 1. counts, shares, medians, ladder -------------------------------
nD = 4;  nA = 4;  rungs = [3 2 1];
OK = false(nD, nA, 3);
OK(:,:,1) = true;  OK(4,4,1) = false;                  % 15 accepted
OK(:,:,2) = true;  OK(4,4,2) = false;  OK(3,3,2) = false;   % 14
OK(:,1:2,3) = true;                                    % 8
base = [1 2 4];
TF = nan(nD, nA, 3);
for k = 1:3
    for r = 1:nD
        TF(r,:,k) = base(k)*(1 + 0.1*r);
    end
end
WALL = 10*ones(nD, nA, 3);  WALL(:,:,1) = 5;
ATT = reshape(1:16, 4, 4);                              % median 8.5
sheet = fullfile(tmp, 'syn.mat');
write_sheet(sheet, rungs, OK, TF, WALL, ATT, 150, 1710);
D = transfer_difficulty(sheet, op);
P = D.perRung;
ok = chk(ok, isequal(P.nAccepted', [15 14 8]), 'nAccepted = [15 14 8]');
ok = chk(ok, all(P.nCells == 16), 'nCells = 16');
ok = chk(ok, max(abs(P.shareAccepted' - [15 14 8]/16)) < 1e-15, 'shareAccepted = count/16');
ok = chk(ok, abs(P.tfMedianND(1) - 1.2) < 1e-12 && abs(P.tfMinND(1) - 1.1) < 1e-12 ...
             && abs(P.tfMaxND(1) - 1.4) < 1e-12, 'rung 1 TF median 1.2, min 1.1, max 1.4');
ok = chk(ok, abs(P.tfSpread(1) - 0.25) < 1e-12, 'rung 1 spread (1.4 - 1.1)/1.2 = 0.25');
ok = chk(ok, abs(P.wallMedianS(1) - 5) < 1e-12 && abs(P.wallMedianS(2) - 10) < 1e-12, 'wall median');
ok = chk(ok, all(isnan(P.nDirect)) && all(isnan(P.revMedian)), 'no sidecar: nDirect and revMedian NaN');
ok = chk(ok, D.ladder.lowestRungAllCells == 3, 'lowestRungAllCells = 3 at shareAll 0.9');
D85 = transfer_difficulty(sheet, struct('print', false, 'shareAll', 0.85));
ok = chk(ok, D85.ladder.lowestRungAllCells == 2, 'lowestRungAllCells = 2 at shareAll 0.85');
ok = chk(ok, D.ladder.lowestRungAnyCell == 1, 'lowestRungAnyCell = 1');
ok = chk(ok, D.ladder.nFullLadders == 8 && D.ladder.nEmptyCells == 1, 'nFullLadders 8, nEmptyCells 1');
ok = chk(ok, abs(D.ladder.attMedian - 8.5) < 1e-12, 'attMedian = 8.5');

% -------- 2. delta-V against SI rocket equation ----------------------------
OKc = true(2, 2, 1);
tfND = 0.5;
sheet2 = fullfile(tmp, 'dv.mat');
write_sheet(sheet2, 5, OKc, tfND*ones(2,2), ones(2,2), ones(2,2), 150, 1710);
Dv = transfer_difficulty(sheet2, op);
tSec = tfND*382981.289129055;
g0 = 9.80665;
mdot = 5/(1710*g0);
mf = 150 - mdot*tSec;
dvRef = 1710*g0/1000*log(150/mf);
ok = chk(ok, abs(Dv.perRung.dVmedianKms - dvRef) < 1e-9*dvRef, ...
         sprintf('dV %.9f km/s vs SI rocket equation %.9f', Dv.perRung.dVmedianKms, dvRef));

% -------- 3. roughness -----------------------------------------------------
ok = chk(ok, all(P.roughness(1:2) >= 0), 'roughness is non-negative');
Dc = transfer_difficulty(sheet2, op);
ok = chk(ok, Dc.perRung.roughness == 0, 'constant sheet: roughness 0 (one-rung sheet handled)');
TFs = ones(3, 3);  TFs(1,1) = 2;
sheet3 = fullfile(tmp, 'spike.mat');
write_sheet(sheet3, 5, true(3,3,1), TFs, ones(3,3), ones(3,3), 150, 1710);
Ds = transfer_difficulty(sheet3, op);
ok = chk(ok, abs(Ds.perRung.roughness - 0.25) < 1e-12, ...
         sprintf('spike sheet roughness = %.6f (hand 0.25)', Ds.perRung.roughness));

% 3b. distinct-neighbour rule on degenerate grids
T14 = [1 2 3 4];
sh = fullfile(tmp, 'g14.mat');
write_sheet(sh, 5, true(1,4,1), T14, ones(1,4), ones(1,4), 150, 1710);
Dg = transfer_difficulty(sh, op);
ok = chk(ok, abs(Dg.perRung.roughness - 0.25) < 1e-12, ...
         sprintf('1 x 4 sheet roughness %.6f (hand 0.25)', Dg.perRung.roughness));
sh = fullfile(tmp, 'g22.mat');
write_sheet(sh, 5, true(2,2,1), [1 2; 3 4], ones(2,2), ones(2,2), 150, 1710);
Dg = transfer_difficulty(sh, op);
ok = chk(ok, abs(Dg.perRung.roughness - 0.3125) < 1e-12, ...
         sprintf('2 x 2 sheet roughness %.6f (hand 0.3125)', Dg.perRung.roughness));

% 3c. >= 2 accepted neighbours: row ends (1 neighbour) and an isolated cell are excluded
OKi = false(5,5,1);  OKi(1,2:4) = true;  OKi(4,4) = true;
Ti = ones(5,5);  Ti(1,2) = 1;  Ti(1,3) = 2;  Ti(1,4) = 4;  Ti(4,4) = 100;
sh = fullfile(tmp, 'iso.mat');
write_sheet(sh, 5, OKi, Ti, ones(5,5), ones(5,5), 150, 1710);
Di = transfer_difficulty(sh, op);
ok = chk(ok, abs(Di.perRung.roughness - 0.25) < 1e-12, ...
         sprintf('only the middle cell has 2 neighbours: roughness %.6f (hand 0.25)', Di.perRung.roughness));

% 3d. unsorted rungs: the smallest thrust, not the last index
OKu = false(2,2,3);  OKu(:,:,1) = true;  OKu(:,:,2) = true;  OKu(1,1,3) = true;
sh = fullfile(tmp, 'uns.mat');
write_sheet(sh, [2 5 1], OKu, ones(2,2,3), ones(2,2,3), ones(2,2), 150, 1710);
Du = transfer_difficulty(sh, op);
ok = chk(ok, Du.ladder.lowestRungAllCells == 2 && Du.ladder.lowestRungAnyCell == 1, ...
         'rungs [2 5 1]: lowestRungAllCells 2, lowestRungAnyCell 1 (smallest thrust)');

% 3e. attempts: ATT = 0 means never attempted
sh = fullfile(tmp, 'att.mat');
write_sheet(sh, 5, true(2,2,1), ones(2,2), ones(2,2), [0 1; 3 5], 150, 1710);
Da = transfer_difficulty(sh, op);
ok = chk(ok, Da.ladder.nAttempted == 3 && Da.ladder.attMedian == 3, ...
         'ATT [0 1 3 5]: nAttempted 3, attMedian 3 (zeros excluded)');

% -------- 4. revolutions ---------------------------------------------------
nDr = 2;  nAr = 2;  tfD = 3;
nt = 401;  t = linspace(0, tfD, nt);
rad = 0.1;  ang = 2*pi*1.5*t/tfD;
relI = [rad*cos(ang); rad*sin(ang); zeros(1, nt)];
rotM = relI;  rotE = relI;
for kk = 1:nt
    cc = cos(t(kk));  ss = sin(t(kk));
    Rm = [cc ss 0; -ss cc 0; 0 0 1];                    % inertial -> rotating
    rotM(:,kk) = Rm*relI(:,kk) + [1 - muStar; 0; 0];    % arc about the Moon
    rotE(:,kk) = Rm*relI(:,kk) + [-muStar; 0; 0];       % same arc about the Earth
end
DX = zeros(7, nt, nDr, nAr, 1);  DX2 = DX;
for a = 1:nDr
    for b = 1:nAr
        DX(1:3,:,a,b,1)  = rotM;
        DX2(1:3,:,a,b,1) = rotE;
    end
end
DDONE = true(nDr, nAr, 1);  DDONE(2,2) = false;
DTF = tfD*ones(nDr, nAr, 1);
sheet4 = fullfile(tmp, 'rev.mat');
write_sheet(sheet4, 5, true(nDr,nAr,1), ones(nDr,nAr), ones(nDr,nAr), ones(nDr,nAr), 150, 1710);
D0 = transfer_difficulty(sheet4, op);
ok = chk(ok, isnan(D0.perRung.revMedian) && isnan(D0.perRung.nDirect), 'no sidecar file: NaN revs');
DU = zeros(3, nt, nDr, nAr, 1); 
save(fullfile(tmp, 'rev_direct.mat'), 'DX', 'DU', 'DTF', 'DDONE');
Dm = transfer_difficulty(sheet4, struct('print', false, 'primary', 'moon'));
ok = chk(ok, abs(Dm.perRung.revMedian - 1.5) < 1e-3, ...
         sprintf('moon-centred inertial arc: revMedian %.6f (want 1.5)', Dm.perRung.revMedian));
ok = chk(ok, Dm.perRung.nDirect == 3, 'DDONE false cell excluded: nDirect = 3');
De = transfer_difficulty(sheet4, struct('print', false, 'primary', 'earth'));
ok = chk(ok, abs(De.perRung.revMedian - 1.5) > 1e-2, ...
         sprintf('same states counted about the Earth differ: %.4f', De.perRung.revMedian));
DX = DX2; 
save(fullfile(tmp, 'rev_direct.mat'), 'DX', 'DU', 'DTF', 'DDONE');
De2 = transfer_difficulty(sheet4, struct('print', false, 'primary', 'earth'));
ok = chk(ok, abs(De2.perRung.revMedian - 1.5) < 1e-3, ...
         sprintf('earth-centred inertial arc: revMedian %.6f (want 1.5)', De2.perRung.revMedian));
bad = false;
try
    transfer_difficulty(sheet4, struct('print', false, 'primary', 'mars'));
catch
    bad = true;
end
ok = chk(ok, bad, 'unknown primary errors');

% -------- 5. empty sheet ---------------------------------------------------
sheet5 = fullfile(tmp, 'empty.mat');
write_sheet(sheet5, [2 1], false(3,3,2), nan(3,3,2), nan(3,3,2), zeros(3,3), 150, 1710);
D5 = transfer_difficulty(sheet5, op);
ok = chk(ok, all(D5.perRung.nAccepted == 0) && all(isnan(D5.perRung.tfMedianND)) ...
             && all(isnan(D5.perRung.dVmedianKms)) && all(isnan(D5.perRung.roughness)), ...
         'no accepted cell: counts 0, statistics NaN');
ok = chk(ok, isnan(D5.ladder.lowestRungAllCells) && isnan(D5.ladder.lowestRungAnyCell) ...
             && D5.ladder.nFullLadders == 0 && D5.ladder.nEmptyCells == 9, 'no accepted cell: ladder summary');
ok = chk(ok, D5.ladder.nAttempted == 0 && isnan(D5.ladder.attMedian), 'all ATT zero: nAttempted 0, attMedian NaN');

% -------- 6. conjPass ------------------------------------------------------
CP = -ones(nD, nA, 3);
CP(:,:,1) = 1;  CP(1,1,1) = 0;  CP(2,1,1) = 0;          % rung 1: 2 refuted of 15
CP(1,2,1) = -1;                                          % one undecided
CP(:,:,3) = 1;  CP(1,1,3) = 0;                           % rung 3: 1 of 8
Dj = transfer_difficulty(sheet, struct('print', false, 'conjPass', CP));
ok = chk(ok, abs(Dj.perRung.shareRefuted(1) - 2/14) < 1e-12 && Dj.perRung.nUndecided(1) == 1, ...
         'rung 1: shareRefuted 2/14, one undecided');
ok = chk(ok, abs(Dj.perRung.shareRefuted(3) - 1/8) < 1e-12 && Dj.perRung.nUndecided(3) == 0, ...
         'rung 3: shareRefuted 1/8');
ok = chk(ok, isnan(Dj.perRung.shareRefuted(2)) && Dj.perRung.nUndecided(2) == 14, ...
         'rung 2: no decided entry gives NaN, 14 undecided');
ok = chk(ok, ~ismember('shareRefuted', P.Properties.VariableNames), 'no conjPass: no shareRefuted column');

% -------- 7. print ---------------------------------------------------------
evalc('Dp = transfer_difficulty(sheet, struct(''print'', true));');
ok = chk(ok, iscell(Dp.lines) && ~isempty(Dp.lines) && ischar(Dp.lines{1}), 'D.lines is a cell of char');
ok = chk(ok, isempty(evalc('transfer_difficulty(sheet, op);')), 'print false is silent');

if ~ok
    error('test_transfer_difficulty:fail', 'one or more checks failed');
end
fprintf('test_transfer_difficulty: ALL PASS\n');
end

% ------------------------------------------------------------------------
function write_sheet(f, rungs, OK, TF, WALL, ATT, m0kg, ispS)
% WRITE_SHEET  Minimal sheet in the engine's format (fields used here).
szD = size(OK, 1);  szA = size(OK, 2);  nR = numel(rungs);
OK = reshape(OK, szD, szA, nR);  TF = reshape(TF, szD, szA, nR);  WALL = reshape(WALL, szD, szA, nR);
FLYKM = nan(szD, szA, nR);  ACCDZ = FLYKM;  RES = FLYKM;  Z8 = nan(8, szD, szA, nR);
sD = (0:szD-1)/szD;  sA = (0:szA-1)/szA;
meta = struct('muStar', 0.012150585609624, 'lStar', 389703.264829278, 'tStar', 382981.289129055, ...
              'depFamily', 'gto', 'arrFamily', 'dro', 'depParams', struct(), 'arrParams', struct(), ...
              'ispS', ispS, 'm0kg', m0kg, 'N', 400);
save(f, 'rungs', 'sD', 'sA', 'TF', 'FLYKM', 'ACCDZ', 'RES', 'WALL', 'OK', 'Z8', 'ATT', 'meta');
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, msg)
% CHK  Print PASS/FAIL for one check and accumulate.
if cond
    fprintf('  PASS  %s\n', msg);
else
    fprintf('  FAIL  %s\n', msg);
    ok = false;
end
end
