function L = unit_lock(action, f, L)
%% Purpose:
%
%   A PROCESS-HELD lock on a file: java.nio FileChannel.tryLock, which the
%   kernel holds on behalf of this MATLAB process and RELEASES when the
%   process exits or is killed. Measured on this host (2026-09-13): a second
%   MATLAB process is refused while the first holds the lock, acquires it
%   the moment the first releases, and acquires it after the first is
%   killed with SIGKILL.
%
%   This is what ownership of a campaign unit means. The first work queue
%   transferred ownership on heartbeat AGE (a claim untouched for 30 min was
%   "abandoned"), which let a slow-but-alive owner and its reclaimer both
%   compute and both publish the same unit (Astra chain review pass 2, A/B).
%   With a process-held lock a live owner cannot be stolen from: a hung
%   owner is KILLED by the supervisor, the kernel frees the lock, and only
%   then can another worker acquire it.
%
%   Advisory, per process, local filesystem. Not for NFS.
%
%% Inputs:
%
%  action                   char                    'try' | 'release' |
%                                                   'probe'
%  f                        char                    lock file path (created
%                                                   if absent)
%  L                        struct                  for 'release': the
%                                                   handle from 'try'
%
%% Outputs:
%
%  L                        struct                  'try'   -> .held logical,
%                                                   .file, and the Java
%                                                   objects that keep it
%                                                   'probe' -> .held = true
%                                                   if SOME process holds it
%                                                   (this one included)
%                                                   'release' -> []
%
%% Revision History:
%  M. Casey                                                   (c) 09/13/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

switch lower(action)
    case 'try'
        L = struct('held', false, 'file', f, 'raf', [], 'ch', [], 'fl', []);
        raf = java.io.RandomAccessFile(f, 'rw');
        ch = raf.getChannel();
        try
            fl = ch.tryLock();                     % null when another PROCESS holds it
        catch                                      % OverlappingFileLockException: THIS process holds it
            fl = [];
        end
        if isempty(fl)
            ch.close();  raf.close();
            return
        end
        L.held = true;  L.raf = raf;  L.ch = ch;  L.fl = fl;

    case 'release'
        if ~isempty(L) && isstruct(L) && L.held
            try L.fl.release(); catch, end
            try L.ch.close();   catch, end
            try L.raf.close();  catch, end
        end
        L = [];

    case 'probe'
        T = unit_lock('try', f);
        L = struct('held', ~T.held);              % could not take it => someone holds it
        if T.held, unit_lock('release', f, T); end

    otherwise
        error('unit_lock:action', 'unknown action "%s"', action);
end
end
