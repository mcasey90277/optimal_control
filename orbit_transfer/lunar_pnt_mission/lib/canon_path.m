function c = canon_path(p)
%% Purpose:
%
%   Canonical absolute path of a file or folder that need not exist: a
%   relative path is taken relative to MATLAB's pwd (NOT the JVM working
%   directory, which java.io.File would use and which does not follow cd),
%   then '.', '..' and symbolic links are resolved (java getCanonicalPath).
%   Shared by export_leg1_transfers (output-folder fence), leg1_movie and
%   full_mission_movie (protected-file fence).
%
%   Self-demo (nargin == 0): a relative and a '..' path.
%
%% Inputs:
%
%  p                        char                    path
%
%% Outputs:
%
%  c                        char                    canonical absolute path
%
%% Revision History:
%  M. Casey  final fix wave: moved out of export_leg1_transfers (c) 11/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    fprintf('canon_path self-demo (pwd %s):\n  x.mp4 -> %s\n  a/../b.mp4 -> %s\n', pwd, ...
            canon_path('x.mp4'), canon_path(fullfile('a', '..', 'b.mp4')));
    return
end
p = char(p);
if ~java.io.File(p).isAbsolute(), p = fullfile(pwd, p); end
c = char(java.io.File(p).getCanonicalPath());
end
