function h = file_sha256(f)
%% Purpose:
%
%   Hex SHA-256 of a file's BYTES (the recipe of leg1_inputs' package
%   identity), computed without loading the file as data. Shared by
%   leg1_package_for and combined_timeline.
%
%   Self-demo (nargin == 0): the hash of this file.
%
%% Inputs:
%
%  f                        char                    file path
%
%% Outputs:
%
%  h                        char [1 x 64]           lower-case hex SHA-256
%                                                   (throws file_sha256:
%                                                   noFile when f is not a
%                                                   file)
%
%% Revision History:
%  M. Casey  final fix wave: one hash recipe                  (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    f = [mfilename('fullpath') '.m'];
    h = file_sha256(f);
    fprintf('file_sha256 self-demo: %s  %s\n', h, f);
    return
end
if ~isfile(f)
    error('file_sha256:noFile', 'file not found: %s', f);
end
md  = java.security.MessageDigest.getInstance('SHA-256');
fid = fopen(f, 'r');
if fid < 0
    error('file_sha256:noFile', 'cannot open %s', f);
end
raw = fread(fid, Inf, '*uint8');
fclose(fid);
md.update(raw);
h = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end
