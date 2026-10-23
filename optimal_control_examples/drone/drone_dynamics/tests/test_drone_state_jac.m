function ok = test_drone_state_jac()
%% Purpose:
%
%   drone_state_jac against two independent differentiations of
%   drone_field: a complex step (exact to round-off, possible because the
%   field is complex-step safe) and central finite differences. Also: the
%   generated file exists, and regenerating it reproduces it byte for byte
%   (so nobody has hand-edited the generated body).
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  ok                       logical                 All checks passed
%
%% Revision History:
%  M. Casey                                                   (c) 10/23/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here);
rngWas = rng(11);  restore = onCleanup(@() rng(rngWas));

models = {'kinematic', 'pointmass', 'rigid6dof'};
for mdl = 1:3
    p = drone_params(models{mdl});
    worstCS = 0;  worstFD = 0;
    for trial = 1:10
        x = randn(p.nx, 1);  if p.nx == 12, x(7:8) = 0.4*x(7:8); end
        u = randn(p.nu, 1);
        A = drone_state_jac(x, u, p);
        Acs = zeros(p.nx);  Afd = zeros(p.nx);
        hc = 1e-20;  hd = 1e-6;
        for k = 1:p.nx
            xc = x;  xc(k) = xc(k) + 1i*hc;
            [Fc, Bc] = drone_field(xc, p);
            Acs(:,k) = imag(Fc + Bc*u)/hc;
            xp = x;  xp(k) = xp(k) + hd;  xm = x;  xm(k) = xm(k) - hd;
            [Fp, Bp] = drone_field(xp, p);  [Fm, Bm] = drone_field(xm, p);
            Afd(:,k) = ((Fp + Bp*u) - (Fm + Bm*u))/(2*hd);
        end
        worstCS = max(worstCS, max(max(abs(A - Acs))));
        worstFD = max(worstFD, max(max(abs(A - Afd))));
    end
    ok = chk(ok, isequal(size(A), [p.nx p.nx]), sprintf('%s: A is %dx%d', p.model, p.nx, p.nx));
    ok = chk(ok, worstCS < 1e-12, sprintf('%s: A vs complex step (%.1e)', p.model, worstCS));
    ok = chk(ok, worstFD < 1e-6, sprintf('%s: A vs central differences (%.1e)', p.model, worstFD));
end

genFile = fullfile(here, 'drone_state_jac_rigid6dof_gen.m');
ok = chk(ok, isfile(genFile), 'the generated rigid6dof Jacobian file exists');
before = fileread(genFile);
restoreGen = onCleanup(@() writeBytes(genFile, before));   % original bytes return on ANY exit path
gen_drone_jacs();
after = fileread(genFile);
ok = chk(ok, strcmp(before, after), 'regenerating reproduces the committed file byte for byte');

if ok, fprintf('TEST_DRONE_STATE_JAC: ALL PASS\n'); else, fprintf('TEST_DRONE_STATE_JAC: FAILURE\n'); end
end

% ------------------------------------------------------------------------
function ok = chk(ok, cond, label)
%% Purpose:
%
%   Accumulate one labelled pass/fail.
%
if cond, fprintf('  PASS  %s\n', label);
else,    fprintf('  FAIL  %s\n', label);  ok = false;
end
end

% ------------------------------------------------------------------------
function writeBytes(fileName, txt)
%% Purpose:
%
%   Write the text back to the file exactly as given (used to restore the
%   committed generated file after the regeneration check).
%
fid = fopen(fileName, 'w');
fwrite(fid, txt, 'char');
fclose(fid);
end
