function bak = guard_catalog_overwrite(catFile, sweepOn)
%% Purpose:
%
%   The check build_70mN_library runs before its PACKAGE stage rewrites a
%   catalog file. The second-order sweep's measurements (conj_interior,
%   conj_near_miss, h6_margin, lift_margin, ...) live ONLY in the catalog:
%   catalog .mat files are gitignored, and second_order_pass's .bak_2nd
%   backup is taken before its writeback, not after. So packaging over a
%   catalog that carries them, with the sweep stage off, destroys hours of
%   measurement with no way back -- which the chain script's shipped
%   switches (package on, sweep off) would have done on a plain run
%   (found 2026-09-11).
%
%   Refuses in that case with a named error. Otherwise copies the existing
%   catalog to a time-stamped backup beside it and returns its path, so
%   every overwrite is reversible.
%
%% Inputs:
%
%  catFile                  char                    the catalog .mat the
%                                                   package stage will write
%  sweepOn                  logical                 true when the chain's
%                                                   sweep stage will rerun
%                                                   and rewrite the fields
%
%% Outputs:
%
%  bak                      char                    backup path, '' when no
%                                                   catalog existed
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

bak = '';
if ~isfile(catFile), return, end
L = load(catFile);
fn = fieldnames(L);
hasSecond = any(cellfun(@(f) isstruct(L.(f)) && isfield(L.(f), 'second_order'), fn));
if hasSecond && ~sweepOn
    error('guard_catalog_overwrite:wouldStripSecondOrder', ...
          ['%s carries the second-order sweep''s measurements, and packaging would erase them ' ...
           'with the sweep stage off. Turn run.sweep on, or package into a different outDir'], catFile);
end
bak = sprintf('%s.bak_%s', catFile, datestr(now, 'yyyymmdd_HHMMSSFFF'));
copyfile(catFile, bak);
end
