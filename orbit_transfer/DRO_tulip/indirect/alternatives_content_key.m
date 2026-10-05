function key = alternatives_content_key(c)
%% Purpose:
%
%   A fingerprint of a catalog's ALTERNATIVES TABLE: for every row, in table
%   order, its departure and arrival phases, its optimality status and its
%   eight numbers. The status audit stores this key next to
%   catalog_content_key, so an audit of the alternatives cannot be credited
%   to a table with other rows, another order or another status.
%
%  ASSUMPTIONS / NOTES:
%
% • Same hashing as catalog_content_key (java MD5 over the raw double bytes).
% • Labels (status_reason, source, gate margins) are left out on purpose:
%   the audit re-derives the status, it does not vouch for the prose.
% • Bit-exact: a costate changed in its last bit is another table.
% • The BORDERLINE flag (relabel_borderline) is bound: it changes what the
%   audit accepts (a lower-bound status). A false flag hashes as an absent
%   one, so the keys of tables without the field are unchanged.
%
%% Inputs:
%
%  c                        struct                  catalog; .alternatives
%                                                   [1 x n] rows with .sD .sA
%                                                   .status .z8 [8 x 1]
%
%% Outputs:
%
%  key                      char [1 x 32]           MD5 of the table, hex;
%                                                   '' when there is none
%
%% Revision History:
%  M. Casey                                                   (c) 10/04/2026
%  M. Casey  borderline flag bound (true rows only)                10/05/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

key = '';
if ~isfield(c, 'alternatives') || isempty(c.alternatives), return, end
A = c.alternatives;
md = java.security.MessageDigest.getInstance('MD5');
md.update(typecast(double(numel(A)), 'uint8'));
for k = 1:numel(A)
    md.update(typecast(double([A(k).sD; A(k).sA; A(k).status; A(k).z8(:)]), 'uint8'));
    if isfield(A, 'borderline') && isequal(A(k).borderline, true)
        md.update(uint8('borderline'));      % a lower-bound status; false hashes as absent (old keys unchanged)
    end
end
key = lower(reshape(dec2hex(typecast(md.digest(), 'uint8'), 2).', 1, []));
end
