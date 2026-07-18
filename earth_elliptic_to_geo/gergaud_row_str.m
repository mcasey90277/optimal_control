function s = gergaud_row_str(row)
% GERGAUD_ROW_STR  Format a gergaud_row() row struct as a fixed-width text
% block matching the Haberkorn-Martinon-Gergaud (JGCD 27(6), 2004) Table-3
% layout.
%
% Pure function: no file I/O, no side effects, deterministic given row.
% When row.certified is false, an UNCERTIFIED banner line (carrying
% row.note) is PREPENDED, so a printed row can never be mistaken for a
% certified result even if only skimmed. Some CERTIFIED rows also carry a
% legitimate caveat in row.note (e.g. the 0.5 N anchor-free R0-law tfmin
% estimate, a custom-endpoint research-probe scope note, or the I-1
% PSR-skipped-for-custom-endpoints note) -- for those, a "NOTE: <note>"
% line is PREPENDED instead, so the caveat is never silently dropped just
% because the row happens to be certified. The two banners are mutually
% exclusive (a row is either uncertified-with-banner or certified-with-
% NOTE-if-any), so a note is never printed twice.
%
% INPUTS:
%   row - struct produced by gergaud_row(), with fields:
%     .thrustN .ctf .tf_ND .tfmin_h .m_f_kg .prop_kg .dV_kms .switches
%     .revs .revs_paper .edge .incl_deg .defect .certified .note        [scalars/char]
%
% OUTPUTS:
%   s - char row: header line + fixed-width data line + defect line,
%       with an "UNCERTIFIED -- <note>" banner prepended when
%       ~row.certified, OR a "NOTE: <note>" line prepended when certified
%       but row.note is non-empty                                       [1xN char]
%
% REFERENCES:
%   [1] Haberkorn, Martinon, Gergaud, "Low Thrust Minimum-Fuel Orbital
%       Transfer: A Homotopic Approach," JGCD 27(6), 2004, Table 3.

if nargin < 1 || ~isstruct(row)
    error('gergaud_row_str:badInput', 'row must be a struct (see gergaud_row)');
end

hdr = sprintf('%-8s %8s %8s %10s %10s %10s %6s %8s %8s %8s %8s %10s\n', ...
    'T[N]', 'tf/tfm', 'tf[ND]', 'tfmin[h]', 'm_f[kg]', 'prop[kg]', 'dV', ...
    'sw', 'revs', 'rev_pap', 'edge', 'incl[deg]');

body = sprintf('%-8.4g %8.4f %8.3f %10.4f %10.2f %10.2f %6.3f %8d %8.3f %8.3f %8.3f %10.4f\n', ...
    row.thrustN, row.ctf, row.tf_ND, row.tfmin_h, row.m_f_kg, row.prop_kg, ...
    row.dV_kms, row.switches, row.revs, row.revs_paper, row.edge, row.incl_deg);

defline = sprintf('defect=%.3e\n', row.defect);

s = [hdr, body, defline];

if ~row.certified
    banner = sprintf('UNCERTIFIED \x2014 %s\n', row.note);
    s = [banner, s];
elseif isfield(row, 'note') && ~isempty(row.note)
    % Certified but footnoted (e.g. anchor-free tfmin estimate, custom-
    % endpoint research-probe scope, or an I-1 PSR-skipped-for-custom-
    % endpoints note): print the caveat so a certified row is never
    % mistaken for a fully-clean certified result (DESIGN Sec 6).
    noteLine = sprintf('NOTE: %s\n', row.note);
    s = [noteLine, s];
end

end
