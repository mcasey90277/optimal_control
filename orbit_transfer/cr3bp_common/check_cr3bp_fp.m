function check_cr3bp_fp(Scached, fpNow, file, tag)
% CHECK_CR3BP_FP  Fail-loud cache-fingerprint guard (earth-campaign pattern).
%
% Compares a loaded cache struct's .fp against the current fingerprint:
%   - no .fp at all       -> WARN (legacy cache, trusted under matching tag)
%   - field only in fpNow -> WARN (schema evolution, compatible)
%   - field on both sides with different values -> ERROR naming field + file
%
% INPUTS:
%   Scached - struct loaded from the cache file [struct]
%   fpNow   - current fingerprint (cr3bp_fingerprint) [struct]
%   file    - cache path, for messages [char]
%   tag     - run tag, for messages [char]
% OUTPUTS: (none) - warns or errors
% REFERENCES:
%   [1] earth_elliptic_to_geo/direct/core/homotopy_mee.m>check_cache_fp (the
%       precedent); [2] spec 2026-07-21-ladder-prep-design.md sec 2.
if ~isfield(Scached, 'fp')
    warning('check_cr3bp_fp:noFingerprint', ...
        '%s has no config fingerprint (legacy cache) -- trusting under tag ''%s''', file, tag);
    return;
end
fn = fieldnames(fpNow);
for k = 1:numel(fn)
    f = fn{k};
    if ~isfield(Scached.fp, f)
        warning('check_cr3bp_fp:schemaOlder', ...
            '%s: fingerprint field ''%s'' absent from cache (schema evolution) -- trusting', file, f);
        continue;
    end
    if ~isequal(Scached.fp.(f), fpNow.(f))
        error('check_cr3bp_fp:mismatch', ...
            ['fingerprint mismatch in %s: field ''%s'' differs from the current ' ...
             'config -- stale/foreign cache under tag ''%s''; delete it or use a new tag'], ...
            file, f, tag);
    end
end
end
