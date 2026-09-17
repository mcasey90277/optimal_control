function out = unit_lock(action, f, token)
%% Purpose:
%
%   A PROCESS-HELD lock on a file, as a lifecycle-controlled capability:
%   java.nio FileChannel.tryLock (POSIX fcntl record locking on macOS),
%   which the kernel holds for this MATLAB process until it releases the
%   lock or dies. Measured on this host (2026-09-13): a second MATLAB
%   process is refused while the first holds it, acquires it the moment
%   the first releases, and acquires it after the first is killed -9.
%
%   THE REGISTRY IS THE AUTHORITY, NOT THE CALLER'S STRUCT. 'try' hands
%   back a handle {file, key, token}; every later question ('holds',
%   'release', 'beat' in the queue) is answered by a process-wide registry
%   that maps the file's IDENTITY (device:inode) to the live Java objects
%   and the token of the CURRENT holder. So a released handle answers
%   "not held" no matter what its copy says, a repeated release is a
%   no-op, and a stale handle cannot act on a newer holder's lock. Astra's
%   pass-3 reproduction -- release, let another process claim, publish
%   with the old struct -- is closed by construction.
%
%   ONE CHANNEL PER FILE PER PROCESS, ENFORCED. POSIX releases ALL of a
%   process's locks on a file when ANY descriptor on that file is closed,
%   and Java documents the same for its channels. So: the registry is
%   consulted BEFORE any open, keyed by identity (not the path string, so
%   aliases cannot slip past); the registry is mlock'ed so `clear
%   functions` / `clear all` cannot wipe it while locks live; and if Java
%   ever reports an overlapping lock (the invariant already broken), that
%   second channel is PARKED, never closed, and the call errors -- closing
%   it is exactly what would drop the live lock (measured, pass 3).
%
%   Advisory, per process, local filesystem. Not for NFS. Never delete or
%   replace a lock file during a campaign: a lock follows the inode.
%
%% Inputs:
%
%  action                   char                    'try' | 'holds' |
%                                                   'release' | 'probe' |
%                                                   'count'
%  f                        char                    lock file path (created
%                                                   if absent)
%  token                    char                    for 'holds'/'release':
%                                                   the token from 'try'
%
%% Outputs:
%
%  out                      'try'     -> struct .held .file .key .token
%                           'holds'   -> logical: this process holds f
%                                        under this token, and the Java
%                                        lock is still valid
%                           'release' -> logical: it was ours and is now
%                                        released (false = no-op)
%                           'probe'   -> struct .held (some process holds
%                                        it) .mine
%                           'count'   -> number of locks this process holds
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

persistent R parked
if isempty(R)
    R = containers.Map('KeyType', 'char', 'ValueType', 'any');
    parked = {};
end

switch lower(action)
    case 'try'
        out = struct('held', false, 'file', f, 'key', '', 'token', '');
        if ~isfile(f), fid = fopen(f, 'a'); if fid >= 0, fclose(fid); end, end
        key = identity(f);  out.key = key;
        if isKey(R, key), return, end            % we hold it: refuse WITHOUT opening
        raf = java.io.RandomAccessFile(f, 'rw');
        ch = raf.getChannel();
        try
            fl = ch.tryLock();                   % null when another PROCESS holds it
        catch ME
            % OverlappingFileLockException: this JVM already holds it and
            % the registry did not know -- the invariant is broken. Park
            % the channel (closing it would free the live lock) and say so.
            parked{end+1} = struct('raf', raf, 'ch', ch);
            error('unit_lock:overlap', ...
                  ['this process already holds a lock on %s that the registry does not know about ' ...
                   '(%s); the second channel has been parked, not closed'], f, ME.message);
        end
        if isempty(fl)
            ch.close();  raf.close();           % we hold nothing on this file: safe to close
            return
        end
        token = char(java.util.UUID.randomUUID());
        R(key) = struct('token', token, 'file', f, 'raf', raf, 'ch', ch, 'fl', fl);
        mlock                                    % the registry must outlive any `clear`
        out.held = true;  out.token = token;

    case 'holds'
        % BY TOKEN: a lock file that was deleted or replaced under us still
        % has a live entry, and must still answer (and be releasable)
        [key, e] = findToken(R, token);
        out = ~isempty(key) && e.fl.isValid();

    case 'release'
        out = false;
        [key, e] = findToken(R, token);
        if isempty(key), return, end                % not ours (or already released)
        remove(R, key);
        try e.fl.release(); catch, end
        try e.ch.close();   catch, end
        try e.raf.close();  catch, end
        if R.Count == 0 && isempty(parked), munlock; end
        out = true;

    case 'probe'
        if isfile(f) && isKey(R, identity(f)), out = struct('held', true, 'mine', true);  return, end
        T = unit_lock('try', f);
        out = struct('held', ~T.held, 'mine', false);   % could not take it => another process holds it
        if T.held, unit_lock('release', f, T.token); end

    case 'count'
        out = R.Count;

    otherwise
        error('unit_lock:action', 'unknown action "%s"', action);
end
end

% ------------------------------------------------------------------------
function [key, e] = findToken(R, token)
% FINDTOKEN  The registry entry holding this token ('' / [] if none). Tokens
% are UUIDs, so the match is exact.  INPUTS: R map; token.  OUTPUTS: key; e.
key = '';  e = [];
if ~ischar(token) || isempty(token), return, end
ks = keys(R);
for k = 1:numel(ks)
    if strcmp(R(ks{k}).token, token), key = ks{k};  e = R(ks{k});  return, end
end
end

% ------------------------------------------------------------------------
function key = identity(f)
% IDENTITY  device:inode of an existing file, so every alias of it (relative
% path, symlink, hard link, case variant on APFS) maps to one registry key.
% Falls back to the canonical path if the attributes are unavailable.
% INPUTS: f.  OUTPUTS: key char.
p = java.io.File(f).toPath();
try
    ino = java.nio.file.Files.getAttribute(p, 'unix:ino');
    dev = java.nio.file.Files.getAttribute(p, 'unix:dev');
    key = sprintf('%d:%d', double(dev), double(ino));
catch
    key = char(java.io.File(f).getCanonicalPath());
end
end
