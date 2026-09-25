function S = video_scenes()
%% Purpose:
%
%   The shot list of the conjugate-point video: one entry per narration
%   scene (narration/sceneNN.txt), each a list of BEATS. A beat fires when
%   the narrator reaches its anchor phrase ('' = the start of the scene);
%   the director converts the phrase's position in the text into a time in
%   that scene's audio, so the pictures follow whatever voice is used.
%
%   Beat operations (applied in order):
%     {'preset', k}              load explorer preset k
%     {'extremal', k}            select extremal k (skipped if absent)
%     {'delta', d}               slope change
%     {'scaled', tf}             difference strip divided by delta or not
%     {'problem', {F,a,b,ya,yb,pRange}}   type a problem and solve
%     {'view', name}             'full' or a panel: curves difference
%                                jacobi shooting mode deltaJ
%     {'caption', txt}           caption bar (TeX); '' removes it
%     {'card', {title, sub}}     full-screen title card
%     {'shrink'}                 animate delta -> 0 across the beat
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  S                        struct array            .file (narration txt)
%                                                   .beats (struct array:
%                                                   .anchor .ops)
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

cat4 = {'y*sqrt(1+yp^2)'};
S = struct('file', {}, 'beats', {});

% 1. cold open: the soap film
S(end+1) = scene('scene01.txt', { ...
  '',                           {{'preset', 4}, {'extremal', 1}, {'view', 'full'}, ...
                                 {'caption', 'A soap film between two rings takes the least area'}}
  'you get two answers',        {{'extremal', 2}, {'view', 'curves'}, ...
                                 {'caption', 'Two curves, one equation'}}
  'Nature only ever builds',    {{'extremal', 1}, {'caption', 'Nature only builds one of them'}}
  'Today we''ll see why',       {{'card', {'Conjugate Points', ...
                                  'when the Euler-Lagrange equation is not enough'}}} });

% 2. setup: extremals are only candidates
S(end+1) = scene('scene02.txt', { ...
  '',                                   {{'preset', 1}, {'view', 'full'}, ...
                                         {'caption', 'Minimize  J[y] = \int_a^b F(t, y, y'') dt'}}
  'The Euler-Lagrange equation finds',  {{'caption', 'Euler-Lagrange equation  \rightarrow  extremals (candidates)'}}
  'But an extremal is like',            {{'caption', 'Like f''(x) = 0:  minimum, maximum, or saddle?'}}
  'For curves, the second derivative',  {{'caption', 'For curves: do neighbouring extremals come back?'}} });

% 3. neighbours that come back
S(end+1) = scene('scene03.txt', { ...
  '',                                       {{'preset', 1}, {'view', 'curves'}, ...
                                             {'caption', 'F = y''^2 - y^2   on  [0, 3]'}}
  'The orange curves are its neighbors',    {{'caption', 'Blue: extremal      Orange: neighbours (same start, new slope)'}}
  'Now stretch the interval out to four',   {{'preset', 2}, {'view', 'curves'}, {'caption', 'Now  b = 4'}}
  'right here, at t equals pi',             {{'caption', 'They all cross again at  t = \pi'}}
  'That point, where the neighboring',      {{'caption', 'Conjugate point: where neighbouring extremals refocus'}}
  'The strip underneath',                   {{'view', 'difference'}, ...
                                             {'caption', 'Neighbour minus extremal:  zeros = crossings'}}
  'And if I divide each difference',        {{'scaled', true}, ...
                                             {'caption', 'Divided by \delta:  every curve lands on the Jacobi field h(t)'}}
  'Its first zero is the conjugate point',  {{'view', 'jacobi'}, ...
                                             {'caption', 'First zero of  h(t) = sin t  is the conjugate point,  \pi'}} });

% 4. why it matters
S(end+1) = scene('scene04.txt', { ...
  '',                                       {{'scaled', false}, {'preset', 2}, {'view', 'full'}, ...
                                             {'caption', 'Jacobi: a conjugate point inside (a, b)  \Rightarrow  not a minimum'}}
  'This bottom panel bends the extremal',   {{'view', 'deltaJ'}, {'caption', 'Cost change along the worst direction'}}
  'the curve points down',                  {{'caption', 'b = 4:  \DeltaJ < 0,  a cheaper nearby path exists'}}
  'Back at three',                          {{'preset', 1}, {'view', 'deltaJ'}, ...
                                             {'caption', 'b = 3:  \DeltaJ > 0,  a genuine local minimum'}}
  'The app counts this a second',           {{'preset', 2}, {'view', 'mode'}, ...
                                             {'caption', 'Morse:  # cost-lowering directions  =  # conjugate points'}} });

% 5. the definition is a limit
S(end+1) = scene('scene05.txt', { ...
  '',                                   {{'preset', 5}, {'extremal', 2}, {'view', 'curves'}, ...
                                         {'caption', 'Pendulum:  F = y''^2/2 + cos y'}}
  'With a big nudge',                   {{'caption', 'Big nudges cross NEAR the conjugate point, not on it'}}
  'Watch what happens as I shrink',     {{'caption', 'Shrinking  \delta \rightarrow 0'}, {'shrink'}}
  'So a conjugate point really lives',  {{'caption', 't_c = 3.301 is a limit: infinitely close neighbours'}}
  'And the rescaled differences',       {{'delta', 0.6}, {'scaled', true}, {'view', 'difference'}, ...
                                         {'caption', '(y_\delta - y_0) / \delta  \rightarrow  Jacobi field h(t)'}} });

% 6. back to the soap film
S(end+1) = scene('scene06.txt', { ...
  '',                                   {{'scaled', false}, {'preset', 4}, {'extremal', 1}, {'view', 'curves'}, ...
                                         {'caption', 'Shallow curve: no conjugate point'}}
  'the cost goes up in every direction', {{'view', 'deltaJ'}, {'caption', 'Shallow:  \DeltaJ > 0,  a minimum'}}
  'The deep curve',                     {{'extremal', 2}, {'view', 'curves'}, ...
                                         {'caption', 'Deep curve: neighbours refocus at  t = 0.153'}}
  'Bend it, and the area drops',        {{'view', 'deltaJ'}, {'caption', 'Deep:  \DeltaJ < 0,  a saddle'}}
  'Now pull the rings apart',           {{'problem', [cat4, {-0.6, 0.6, 1, 1, [-8 2]}]}, {'extremal', 2}, ...
                                         {'view', 'curves'}, {'caption', 'Half-gap  0.60'}}
  'In the shooting panel',              {{'view', 'shooting'}, ...
                                         {'caption', 'Shooting function: two roots = two catenaries'}}
  'At a half gap of about',             {{'problem', [cat4, {-0.66, 0.66, 1, 1, [-8 2]}]}, {'view', 'shooting'}, ...
                                         {'caption', 'Half-gap  0.66:  the two solutions nearly merge'}}
  'Go any wider',                       {{'problem', [cat4, {-0.67, 0.67, 1, 1, [-8 2]}]}, {'view', 'shooting'}, ...
                                         {'caption', 'Half-gap  0.67:  no solution at all'}}
  'Physically, that''s the moment',     {{'caption', 'The soap film pops'}} });

% 7. wrap-up
S(end+1) = scene('scene07.txt', { ...
  '',                                       {{'preset', 4}, {'extremal', 2}, {'view', 'full'}, ...
                                             {'caption', 'Four views of one idea'}}
  'Neighbors crossing',                     {{'view', 'curves'}, {'caption', '1.  Neighbours crossing'}}
  'The Jacobi field hitting zero',          {{'view', 'jacobi'}, {'caption', '2.  The Jacobi field hitting zero'}}
  'The cost actually going down',           {{'view', 'deltaJ'}, {'caption', '3.  The cost going down'}}
  'And the number of those directions',     {{'view', 'mode'}, {'caption', '4.  The Morse count'}}
  'It''s the same test engineers use',      {{'card', {'The same test, bigger problems', ...
                                              'certifying optimal spacecraft trajectories'}}}
  'Now continue your night walk',           {{'view', 'full'}, {'caption', ''}} });
end

% ---------------------------------------------------------------------------
function s = scene(file, rows)
% SCENE  One scene record. INPUTS: narration file name; rows {anchor, ops}.
% OUTPUTS: s struct with .file and .beats.
beats = struct('anchor', rows(:,1), 'ops', rows(:,2));
s = struct('file', file, 'beats', beats);
end
