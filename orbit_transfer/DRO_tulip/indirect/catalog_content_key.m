function key = catalog_content_key(c)
%% Purpose:
%
%   A fingerprint of WHAT a phase catalog holds: which cells have a root, the
%   flight time of each and its eight numbers. An audit stores the key of the
%   catalog it read; whoever later relies on that audit checks the key, so an
%   audit cannot be credited to a catalog with the same number of entries
%   and other contents (Astra review 2026-09-19).
%
%  ASSUMPTIONS / NOTES:
%
% • Labels are left out on purpose (family_index, notes, verdict fields): they
%   are rewritten after the audit and do not change what was certified.
% • Only cells with a solution count; what sits elsewhere in tf_nd does not.
% • Bit-exact: a costate changed in its last bit is another catalog.
%
%% Inputs:
%
%  c                        struct                  phase catalog, c.sheets(k)
%                                                   with .has_solution .tf_nd
%                                                   .z8
%
%% Outputs:
%
%  key                      char [1 x 32]           MD5 of the content, hex
%
%% Revision History:
%  M. Casey                                                   (c) 09/20/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

md = java.security.MessageDigest.getInstance('MD5');
for ks = 1:numel(c.sheets)
    s = c.sheets(ks);
    has = logical(s.has_solution);
    md.update(typecast(double(size(has)), 'uint8'));
    md.update(uint8(has(:)));
    md.update(typecast(double(s.tf_nd(has)), 'uint8'));
    md.update(typecast(double(s.z8(:)), 'uint8'));
end
key = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end
