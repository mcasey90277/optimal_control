% TEST_SOLVER_FORK_PARITY  Guard the tulip/ELFO solver fork against silent drift.
%
% WHY THIS EXISTS. GTO_tulip's casadi_minfuel_sundman and GTO_ELFO's
% casadi_energy_freetf are a deliberate FORK: ELFO added a 9th state (cScale,
% Betts' slack-state free-time trick), a free horizon and a two-primary clock,
% and that differing quarter is genuinely the substance. Merging them was
% considered on 2026-07-26 and refused -- it would produce a heavily conditional
% solver and put two campaigns' certified results at risk.
%
% But refusing a merge leaves the fork's real hazard unaddressed: a fix landed
% in one copy does not reach the other, and NOTHING NOTICES. That is exactly the
% failure mode that cost this repo a silently-dead feature for eleven days when
% PSR/lib's vendored copies drifted from upstream. Dissolving PSR/lib removed
% that hazard; this fork carries the same one, so it gets the same treatment.
%
% WHAT IS GUARDED. Not byte-identity -- these files are supposed to differ.
% Instead: the contiguous runs of code the two still SHARE (>= 4 lines, comments
% and whitespace normalized away). Measured 2026-07-26: 8 runs, 87 lines,
% covering the state box constraints, the warm-start priming, the defect
% assembly, the dual/costate extraction, the output struct and the CasADi
% bootstrap. Each is pinned below by length and CRC.
%
% Note what is NOT in the list any more: the IPOPT options block, formerly the
% single largest shared run at 41 lines, is gone -- it was extracted to
% cr3bp_common/cr3bp_ipopt_opts.m. That is what a resolved duplication looks
% like, and the shrinking of this manifest over time is the goal.
%
% WHEN THIS FAILS, DO NOT JUST RE-PIN. A broken block means one copy changed and
% the other did not. Decide deliberately:
%   * the change belongs in both -> port it, then re-pin;
%   * the change is genuinely ELFO-only (or tulip-only) -> the block has left
%     the shared set; delete its row and say why in the commit.
% Re-pinning without deciding is how a fork silently rots.
%
% INPUTS:  none
% OUTPUTS: none (throws, naming the block that diverged)
%
% REFERENCES:
%   [1] orbit_transfer/CODE_STRUCTURE.md (why the fork was kept).
%   [2] GTO_ELFO/doc/gto_elfo_guide.tex sec 3 (what the 24% buys).
here = fileparts(mfilename('fullpath'));
ot   = fileparts(fileparts(here));                       % orbit_transfer/
fT   = fullfile(ot, 'GTO_tulip', 'direct', 'lib',  'casadi_minfuel_sundman.m');
fE   = fullfile(ot, 'GTO_ELFO',  'direct', 'elfo', 'casadi_energy_freetf.m');
assert(isfile(fT), 'tulip solver not found at %s', fT);
assert(isfile(fE), 'ELFO solver not found at %s',  fE);

% Pinned shared runs: {nLines, crc32, label}. Measured 2026-07-26.
BLOCKS = {
    27, 2515234513, 'opti.set_initial(X, X0);'
    16,  799462778, 'opti.subject_to(X(:) >= lbX(:)); ... state boxes'
    12, 2224613960, '''maxDefect'', max(abs(Dd(:))), ... output struct'
    12, 3055658429, 'lamV = lamDef(4:6, :); ... dual/costate extraction'
     6, 2694191111, 'cpath = getenv(''CASADI_PATH''); ... bootstrap'
     5, 2714178885, 'opti.subject_to(D(:) == 0); ... defect constraints'
     5, 1418932714, 'U = opti.variable(4, nN); ... decision variables'
     4,  291843625, '''nonphysical box '' ... bound-saturation diagnostic'
};

linesT = local_codelines(fT);
linesE = local_codelines(fE);

fprintf('\n=== tulip/ELFO solver fork parity ===\n');
fprintf('  tulip %d code lines, ELFO %d\n', numel(linesT), numel(linesE));

bad = {};
for k = 1:size(BLOCKS,1)
    n = BLOCKS{k,1};  want = BLOCKS{k,2};  lab = BLOCKS{k,3};
    inT = local_hasblock(linesT, n, want);
    inE = local_hasblock(linesE, n, want);
    if inT && inE
        fprintf('  %2d lines  shared    %s\n', n, lab);
    else
        where = 'NEITHER';
        if inT && ~inE, where = 'tulip only -- ELFO copy changed';
        elseif ~inT && inE, where = 'ELFO only -- tulip copy changed';
        end
        fprintf('  %2d lines  DIVERGED  %s   <-- %s\n', n, lab, where);
        bad{end+1} = sprintf('%s  [%s]', lab, where); %#ok<AGROW>
    end
end

if ~isempty(bad)
    error('test_solver_fork_parity:drift', ...
        ['%d shared block(s) no longer match across the tulip/ELFO solver fork:\n' ...
         '  - %s\n\n' ...
         'One copy changed and the other did not. Port the change to both and ' ...
         're-pin, OR delete the row if the divergence is deliberate and say why. ' ...
         'Do NOT re-pin without deciding -- that is how a fork silently rots.'], ...
        numel(bad), strjoin(bad, sprintf('\n  - ')));
end

fprintf('\ntest_solver_fork_parity: ALL PASS (%d shared blocks, %d lines, still in sync)\n', ...
        size(BLOCKS,1), sum([BLOCKS{:,1}]));

% ---------------------------------------------------------------------------
function L = local_codelines(f)
% LOCAL_CODELINES  Code lines with comments and surrounding whitespace removed.
% INPUTS:  f - path to a .m file [char]
% OUTPUTS: L - non-empty normalized code lines [cell of char]
txt = fileread(f);
raw = regexp(txt, '\r?\n', 'split');
L = {};
for k = 1:numel(raw)
    c = raw{k};
    p = strfind(c, '%');            % strip from the first comment character
    if ~isempty(p), c = c(1:p(1)-1); end
    c = strtrim(c);
    if ~isempty(c), L{end+1} = c; end %#ok<AGROW>
end
end

% ---------------------------------------------------------------------------
function tf = local_hasblock(L, n, wantCrc)
% LOCAL_HASBLOCK  True if some n-line window of L has the given CRC32.
% INPUTS:  L - code lines [cell]; n - window length [scalar]; wantCrc [scalar]
% OUTPUTS: tf - logical
tf = false;
for s = 1:numel(L)-n+1
    if local_crc32(strjoin(L(s:s+n-1), sprintf('\n'))) == wantCrc
        tf = true; return
    end
end
end

% ---------------------------------------------------------------------------
function c = local_crc32(str)
% LOCAL_CRC32  CRC-32 (IEEE 802.3) of a char string. Self-contained so the test
% has no toolbox or Java dependency.
% INPUTS:  str - text [char]
% OUTPUTS: c   - CRC-32 [uint32-valued double]
persistent TBL
if isempty(TBL)
    TBL = zeros(1,256,'uint32');
    for b = 0:255
        c0 = uint32(b);
        for bit = 1:8
            if bitand(c0,1), c0 = bitxor(bitshift(c0,-1), uint32(3988292384));
            else,            c0 = bitshift(c0,-1);
            end
        end
        TBL(b+1) = c0;
    end
end
crc = uint32(hex2dec('FFFFFFFF'));
bytes = uint8(unicode2native(str, 'UTF-8'));
for k = 1:numel(bytes)
    idx = double(bitand(bitxor(crc, uint32(bytes(k))), uint32(255))) + 1;
    crc = bitxor(bitshift(crc,-8), TBL(idx));
end
c = double(bitxor(crc, uint32(hex2dec('FFFFFFFF'))));
end
