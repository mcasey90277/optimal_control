function nCopied = adopt_walked_arcs(srcDir, names, arcDir, tag, opts)
%% Purpose:
%
%   Give a phase-torus campaign the arrival arcs that were ALREADY WALKED for
%   its families, so run_phase_torus adopts them instead of walking them
%   again. An arc is 2 to 6 hours of continuation; the ten arcs of the 70 mN
%   library are most of a rebuild's compute, and they are on disk.
%
%   The hand-built campaign kept its arcs in indirect/results/ as
%       arrival_arc_<name>_<dn|up>_long.mat
%   and run_phase_torus looks in <outDir>/arcs for
%       arrival_arc_<tag>_<name>_<dn|up>_long.mat
%   so adoption is a copy under the driver's name, both directions of every
%   named family.
%
%  ASSUMPTIONS / NOTES:
%
% • Nothing is ever overwritten. A destination that already holds the same
%   bytes is left alone (so the call is safe to repeat); one that holds
%   anything else is an error, because that arc belongs to another walk.
% • A missing source is an error by name, never a silent skip: a family
%   without its arcs would quietly drop out of the rebuilt sheet.
%
%% Inputs:
%
%  srcDir                   char                    folder of the walked arcs
%  names                    cellstr                 the anchor (family) names
%  arcDir                   char                    the campaign's arc folder
%                                                   (created if absent)
%  tag                      char                    the campaign's tag
%  opts                     struct (optional)       .print [true]
%
%% Outputs:
%
%  nCopied                  int                     files copied by this call
%
%% Revision History:
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin < 5, opts = struct(); end
doPrint = ~isfield(opts, 'print') || opts.print;
if ~isfolder(arcDir), mkdir(arcDir); end
nCopied = 0;
for k = 1:numel(names)
    for dirn = {'dn', 'up'}
        src = fullfile(srcDir, sprintf('arrival_arc_%s_%s_long.mat', names{k}, dirn{1}));
        dst = fullfile(arcDir, sprintf('arrival_arc_%s_%s_%s_long.mat', tag, names{k}, dirn{1}));
        assert(isfile(src), 'adopt_walked_arcs:missing', ...
               'family %s has no walked %s arc: %s is missing', names{k}, dirn{1}, src);
        if isfile(dst)
            assert(sameBytes(src, dst), 'adopt_walked_arcs:exists', ...
                   '%s already exists and is NOT a copy of %s; it is not overwritten', dst, src);
            continue                                   % already adopted
        end
        copyfile(src, dst);
        nCopied = nCopied + 1;
        if doPrint, fprintf('  adopted %-10s %s  ->  %s\n', names{k}, dirn{1}, dst); end
    end
end
end

% ------------------------------------------------------------------------
function tf = sameBytes(a, b)
% SAMEBYTES  Do two files hold the same bytes?  INPUTS: a; b (paths).
% OUTPUTS: tf.
da = dir(a);  db = dir(b);
tf = da.bytes == db.bytes;
if ~tf, return, end
fa = fopen(a, 'r');  xa = fread(fa, inf, '*uint8');  fclose(fa);
fb = fopen(b, 'r');  xb = fread(fb, inf, '*uint8');  fclose(fb);
tf = isequal(xa, xb);
end
