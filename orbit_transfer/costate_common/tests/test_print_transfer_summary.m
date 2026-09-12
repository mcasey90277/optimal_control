function ok = test_print_transfer_summary()
%% Purpose:
%
%   Tests print_transfer_summary -- the transfer numbers printed ONE way for
%   both consumers: the study script (a fly_transfer flight) and
%   run_dro_tulip (a certify_root certificate), which had been printing the
%   same three numbers in two formats.
%
%   Checks: a flight prints two lines including its clearances; a
%   certificate, which carries no flight, prints one; the numbers are the
%   struct's own; quiet returns the lines without printing; a struct missing
%   a required field is refused by name.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/11/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
addpath(here, fullfile(fileparts(here), 'DRO_tulip', 'indirect'));

flight = struct('tfNd', 4.0151, 'tfDays', 17.7976, 'dvKms', 0.7485, 'propellantKg', 12.2031, ...
                'nSamples', 4197, 'tStar', 382981.289129055, ...
                'admissibility', struct('ok', true, 'reason', 'admissible', ...
                                        'moonKm', 4673.2, 'earthKm', 365000.9));
L = print_transfer_summary(flight, struct('quiet', true));
ok = chk(ok, numel(L) == 2, sprintf('a flight prints two lines (%d)', numel(L)));
ok = chk(ok, contains(L{1}, '0.7485') && contains(L{1}, '12.20') && contains(L{1}, '17.7976'), ...
         'the first line carries t_f, Delta-V and propellant');
ok = chk(ok, contains(L{2}, '4197') && contains(L{2}, '4673') && contains(L{2}, 'admissible'), ...
         'the second carries the samples, the verdict and the clearances');

cert = struct('tfDays', 17.8775, 'dvKms', 0.7512, 'propellantKg', 12.25);
Lc = print_transfer_summary(cert, struct('quiet', true));
ok = chk(ok, numel(Lc) == 1 && contains(Lc{1}, '17.8775'), ...
         'a certificate with no flight prints one line');

out = evalc('print_transfer_summary(flight);');
ok = chk(ok, count(out, newline) == 2, 'printing emits exactly the lines it built');
ok = chk(ok, isempty(evalc('print_transfer_summary(flight, struct(''quiet'', true));')), ...
         'quiet prints nothing');

ok = chk(ok, refuses(@() print_transfer_summary(rmfield(cert, 'dvKms'), struct('quiet', true)), ...
         'print_transfer_summary:fields'), 'a struct with no Delta-V is refused');
ok = chk(ok, refuses(@() print_transfer_summary(struct('dvKms', 1, 'propellantKg', 2), struct('quiet', true)), ...
         'print_transfer_summary:fields'), 'a struct with no time of flight is refused');

if ok, fprintf('TEST_PRINT_TRANSFER_SUMMARY: ALL PASS\n'); else, fprintf('TEST_PRINT_TRANSFER_SUMMARY: FAIL\n'); end
end

function r = refuses(fh, id)
% REFUSES  True when fh throws the named identifier.  INPUTS: fh; id.
% OUTPUTS: r [logical].
try
    fh();  r = false;
catch ME
    r = strcmp(ME.identifier, id);
    if ~r, fprintf('        (threw %s: %s)\n', ME.identifier, ME.message); end
end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
