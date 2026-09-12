function ok = test_movie_phase_sweep()
%% Purpose:
%
%   Tests movie_phase_sweep -- the (departure phase x arrival phase) sweep
%   movie. The point of the function is that everything EXCEPT the
%   trajectory is identical frame to frame, so that is what is tested:
%
%     1. the phase SELECTION refuses what it cannot draw, by name: an
%        off-grid value, both input forms for one axis, a bad order, a bad
%        index, and a selection with no certified pair in it;
%     2. the LOOP ORDER does what it says, and both orders cover the same
%        set of pairs and produce the same colour axis;
%     3. the COLOUR AXIS is [0, max t_f over the selection], not per frame;
%     4. the AXIS BOX is equal-aspect and shared;
%     5. the TITLE's second line is the SAME LENGTH in every frame,
%        including a gap frame -- the constant-width requirement, tested
%        rather than asserted (the titles are returned for this reason);
%     6. a pair with no certified entry still produces a frame;
%     7. .slow holds every frame longer in BOTH outputs, by the same factor,
%        and a non-positive value is refused.
%
%   Rendering is real but small: 2 to 4 frames at about a second each.
%
%% Inputs:  none
%% Outputs: ok [logical]
%
%% Revision History:
%  M. Casey                                                   (c) 09/12/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ok = true;
here = fileparts(fileparts(mfilename('fullpath')));
ind = fullfile(fileparts(here), 'DRO_tulip', 'indirect');
addpath(here, ind);
stem = [tempname '_sweep'];
cleanup = onCleanup(@() cellfun(@(f) delDbl(f), {[stem '.mp4'], [stem '.gif']}));

% ---- 1. selection refuses what it cannot draw ---------------------------
bad = {struct('sD', 0, 'sA', 0.123),        'movie_phase_sweep:offGrid',    'an off-grid arrival value';
       struct('sD', 0, 'sDIdx', 1),         'movie_phase_sweep:bothForms',  'values AND indices for one axis';
       struct('order', 3),                  'movie_phase_sweep:order',      'an order that is neither 1 nor 2';
       struct('sAIdx', 99),                 'movie_phase_sweep:badIndex',   'an index off the grid';
       struct('sDIdx', 2, 'sAIdx', 10),     'movie_phase_sweep:empty',      'a selection with no certified pair';
       struct('sDIdx', 1, 'sAIdx', 2, 'slow', 0), 'movie_phase_sweep:slow',  'a slow factor of zero'};
for k = 1:size(bad, 1)
    o = bad{k,1};  o.outStem = stem;
    try
        movie_phase_sweep(o);
        ok = chk(ok, false, sprintf('%s should have been refused', bad{k,3}));
    catch ME
        ok = chk(ok, strcmp(ME.identifier, bad{k,2}), ...
                 sprintf('%s refused as %s', bad{k,3}, ME.identifier));
    end
end

% ---- 2, 3, 4. orders, colour axis, box -----------------------------------
o1 = movie_phase_sweep(struct('sDIdx', [1 3], 'sAIdx', [2 3], 'order', 1, 'outStem', stem));
o2 = movie_phase_sweep(struct('sDIdx', [1 3], 'sAIdx', [2 3], 'order', 2, 'outStem', stem));
ok = chk(ok, o1.nFrames == 4 && o2.nFrames == 4, 'both orders render every pair (4 frames)');
ok = chk(ok, o1.pairs(1,1) == o1.pairs(2,1) && o1.pairs(1,2) ~= o1.pairs(2,2), ...
         'order 1 holds sD and advances sA');
ok = chk(ok, o2.pairs(1,2) == o2.pairs(2,2) && o2.pairs(1,1) ~= o2.pairs(2,1), ...
         'order 2 holds sA and advances sD');
ok = chk(ok, isequal(sortrows(round(o1.pairs,9)), sortrows(round(o2.pairs,9))), ...
         'both orders cover the same set of pairs');
ok = chk(ok, isequal(o1.clim, o2.clim) && isequal(o1.axLim, o2.axLim), ...
         'colour axis and box do not depend on the order');
ok = chk(ok, o1.clim(1) == 0 && abs(o1.clim(2) - max(o1.tfDays)) < 1e-12, ...
         sprintf('colour axis is [0, max t_f] = [0 %.3f] d, not per frame', o1.clim(2)));
rng_ = diff(o1.axLim, 1, 2);
ok = chk(ok, max(abs(rng_ - rng_(1))) < 1e-12, ...
         sprintf('the shared box is equal-aspect (%.4f ND on every axis)', rng_(1)));

% ---- 5, 6. constant title width, including a gap frame -------------------
% arrival index 10 has no certified entry at any sD in the shipped catalog,
% so this selection mixes a drawn frame with a gap frame
og = movie_phase_sweep(struct('sDIdx', 3, 'sAIdx', [2 10], 'order', 1, 'outStem', stem));
ok = chk(ok, og.nFrames == 2 && og.nCertified == 1 && og.nGaps == 1, ...
         sprintf('a gap still gets a frame: %d certified, %d gap', og.nCertified, og.nGaps));
L2 = cellfun(@(t) numel(t{2}), og.titles);
ok = chk(ok, numel(unique(L2)) == 1, ...
         sprintf('certified and gap titles are the SAME LENGTH (%s chars)', mat2str(unique(L2))));
L1 = cellfun(@(t) numel(t{2}), o1.titles);
ok = chk(ok, numel(unique(L1)) == 1, ...
         sprintf('and every drawn frame''s title is too (%s chars, t_f %.3f .. %.3f d)', ...
                 mat2str(unique(L1)), min(o1.tfDays), max(o1.tfDays)));
ok = chk(ok, contains(og.titles{2}{2}, '-----'), 'the gap title dashes the numbers rather than dropping the fields');

% ---- 7. slow scales both outputs by the same factor ----------------------
oN = movie_phase_sweep(struct('sDIdx', [1 3], 'sAIdx', 2, 'outStem', stem));
oS = movie_phase_sweep(struct('sDIdx', [1 3], 'sAIdx', 2, 'slow', 2, 'outStem', stem));
ok = chk(ok, abs(oS.fps - oN.fps/2) < 1e-12, ...
         sprintf('slow = 2 halves the video rate: %.3g -> %.3g fps', oN.fps, oS.fps));
ok = chk(ok, abs(oS.gifDelay - 2*oN.gifDelay) < 1e-12, ...
         sprintf('and doubles the gif delay: %.4f -> %.4f s', oN.gifDelay, oS.gifDelay));
ok = chk(ok, abs(oS.runSec - 2*oN.runSec) < 1e-12, ...
         sprintf('so the movie runs twice as long: %.2f -> %.2f s for the same %d frames', ...
                 oN.runSec, oS.runSec, oS.nFrames));
ok = chk(ok, oS.nFrames == oN.nFrames, 'slowing holds frames longer, it does not add frames');
vS = VideoReader(oS.mp4);
ok = chk(ok, abs(vS.FrameRate - oS.fps) < 1e-6, ...
         sprintf('the written mp4 carries the slowed rate (%.3g fps)', vS.FrameRate));

% ---- the files exist and are the declared size ---------------------------
v = VideoReader(og.mp4);
ok = chk(ok, v.Height == 720 && v.Width == 1280, ...
         sprintf('frames are exactly 1280x720 (H.264 shears otherwise): %dx%d', v.Width, v.Height));
ok = chk(ok, isfile(og.gif), 'a gif is written beside the mp4');

if ok, fprintf('TEST_MOVIE_PHASE_SWEEP: ALL PASS\n'); else, fprintf('TEST_MOVIE_PHASE_SWEEP: FAIL\n'); end
end

function delDbl(f)
% DELDBL  Delete a file if it exists.  INPUTS: f.  OUTPUTS: none.
if isfile(f), delete(f); end
end

function ok = chk(ok, cond, msg)
% CHK  Print one check.  INPUTS: ok; cond; msg.  OUTPUTS: ok.
if cond, fprintf('  PASS  %s\n', msg); else, fprintf('  FAIL  %s\n', msg); ok = false; end
end
