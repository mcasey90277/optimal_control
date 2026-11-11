function fp = catalog_identity(M)
%% Purpose:
%
%   THE CATALOG IDENTITY of the lunar PNT mission study: the fields of the
%   study's results fingerprint that say which catalog, regime, slot
%   phasing, spacecraft phases and orbit handles a run used. ONE recipe,
%   used by the study (which appends its polish thresholds and polishGates
%   to make the full fingerprint) and by export_mission_transfers (which
%   refuses results whose identity differs from the catalog on the path).
%   MD5s are of a numeric array's bytes with its size prepended (doubles).
%
%   Self-demo (nargin == 0): the identity of the catalog on the path.
%
%% Inputs:
%
%  M                        struct                  mission_inputs() output
%
%% Outputs:
%
%  fp                       struct                  .catalogName .catalogCreated
%                                                   (char) .z8md5 .tfmd5 (char
%                                                   [1 x 32], MD5 of sh.z8 and
%                                                   sh.tf_nd) .thrustN .ispS
%                                                   .m0kg .muStar .lStar .tStar
%                                                   (scalars) .dtau6 [1 x 6]
%                                                   .droPhase0 [1 x 6] .droMd5
%                                                   .tulipMd5 (char [1 x 32],
%                                                   the DRO and tulip handles
%                                                   at 64 phases)
%
%% Revision History:
%  M. Casey  moved out of the study and the export (one recipe) (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    here = fileparts(fileparts(mfilename('fullpath')));
    run(fullfile(here, 'setup_lunar_pnt_paths.m'));
    M = mission_inputs();
end
s64 = (0:63)/64;
fp = struct('catalogName', char(string(M.cat.name)), 'catalogCreated', char(string(M.cat.created)), ...
            'z8md5', md5_of(M.sh.z8), 'tfmd5', md5_of(M.sh.tf_nd), ...
            'thrustN', M.thrustN, 'ispS', M.ispS, 'm0kg', M.m0kg, 'muStar', M.muStar, ...
            'lStar', M.lStar, 'tStar', M.tStar, 'dtau6', M.dtau6, ...
            'droPhase0', M.droPhase0, 'droMd5', md5_of(M.depState(s64)), ...
            'tulipMd5', md5_of(M.arrState(s64)));
if nargin == 0, disp(fp); end
end

% ------------------------------------------------------------------------
function h = md5_of(x)
% MD5_OF  Hex MD5 of a numeric array's bytes (size included).
%   INPUTS: x numeric array.  OUTPUTS: h char [1 x 32].
md = java.security.MessageDigest.getInstance('MD5');
md.update(typecast([double(size(x)), double(x(:)).'], 'uint8'));
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end
