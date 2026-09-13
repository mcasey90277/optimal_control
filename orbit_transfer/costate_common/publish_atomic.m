function publish_atomic(src, dst)
%% Purpose:
%
%   Publish a finished file: move src onto dst in ONE rename(2), replacing
%   any existing dst, so that a reader sees either the old file or the new
%   one and never a partial or missing one. This is the only way a campaign
%   artifact, queue record or heartbeat reaches its final name.
%
%   Why not movefile. Measured on this host (R2026a, APFS, 2026-09-13):
%   MATLAB's movefile onto an existing FILE does rename (the destination
%   takes the source's inode), but that is an observation, not a contract,
%   and movefile onto an existing DIRECTORY moves the source INSIDE it --
%   which is how a claim takeover could have nested one claim in another.
%   java.nio Files.move with ATOMIC_MOVE is DOCUMENTED to rename or throw:
%   it never degrades to a copy, and it refuses a directory destination
%   ("Is a directory", measured). (Astra chain review pass 2, item C.)
%
%% Inputs:
%
%  src                      char                    the finished temporary
%                                                   file, in the SAME folder
%                                                   as dst (same filesystem)
%  dst                      char                    the final name
%
%% Outputs:  none (throws on any failure)
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

assert(isfile(src), 'publish_atomic:src', 'nothing to publish: %s is not a file', src);
assert(~isfolder(dst), 'publish_atomic:dst', 'refusing to publish onto a directory: %s', dst);
[ds, ~] = fileparts(dst);  [ss, ~] = fileparts(src);
assert(strcmp(ds, ss), 'publish_atomic:dir', ...
       'src and dst must share a folder (a rename cannot cross filesystems): %s vs %s', ss, ds);
import java.nio.file.*
try
    Files.move(java.io.File(src).toPath(), java.io.File(dst).toPath(), ...
               [StandardCopyOption.ATOMIC_MOVE, StandardCopyOption.REPLACE_EXISTING]);
catch ME
    error('publish_atomic:move', 'atomic move %s -> %s failed: %s', src, dst, ...
          regexprep(ME.message, '\s+', ' '));
end
assert(isfile(dst) && ~isfile(src), 'publish_atomic:post', 'after the move %s should exist and %s should not', dst, src);
end
