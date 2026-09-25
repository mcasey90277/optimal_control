function P = cov_presets()
%% Purpose:
%
%   The ready-made problems of the conjugate-point explorer and its guide.
%   Each is chosen to show one thing:
%
%     1-3  F = y'^2 - y^2 (harmonic oscillator). Jacobi field sin(t - a),
%          conjugate points at a + k pi. Before pi: a minimiser. Past pi:
%          not. Past 2 pi: two conjugate points, two negative modes.
%     4    F = y sqrt(1 + y'^2) (minimal surface of revolution). Two
%          catenaries join the same ends; the shallow one is a minimiser,
%          the deep one has a conjugate point inside the interval.
%     5    F = y'^2/2 + cos(y) (pendulum action). Nonlinear: a finite slope
%          change crosses NEAR the conjugate point, not at it.
%     6    F = y'^2 + y^2. Q > 0: no conjugate point on any interval.
%     7    F = sqrt(1 + y'^2) (arc length). Straight lines, Jacobi field
%          t - a, never zero: the shortest path.
%
%% Inputs:
%
%  none
%
%% Outputs:
%
%  P                        struct array            .name .F .a .b .ya .yb
%                                                   .pRange .delta .note
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

P = struct('name', {}, 'F', {}, 'a', {}, 'b', {}, 'ya', {}, 'yb', {}, ...
           'pRange', {}, 'delta', {}, 'note', {});
P(end+1) = mk('Oscillator, b = 3 (before pi)', 'yp^2 - y^2', 0, 3, 0, 1, [-12 12], 0.5, ...
              'Jacobi field sin(t): its zero at pi lies beyond b = 3. A minimiser.');
P(end+1) = mk('Oscillator, b = 4 (past pi)', 'yp^2 - y^2', 0, 4, 0, 1, [-5 5], 0.5, ...
              'Every neighbour meets the extremal again at t = pi < b. Not a minimiser.');
P(end+1) = mk('Oscillator, b = 7 (past 2 pi)', 'yp^2 - y^2', 0, 7, 0, 1, [-5 5], 0.5, ...
              'Two conjugate points (pi, 2 pi): two negative modes of the second variation.');
P(end+1) = mk('Minimal surface (two catenaries)', 'y*sqrt(1+yp^2)', -0.5, 0.5, 1, 1, [-8 2], 0.4, ...
              'Two extremals. The shallow one minimises; the deep one has a conjugate point.');
P(end+1) = mk('Pendulum (nonlinear)', 'yp^2/2 + cos(y)', 0, 5, 0, 0.5, [-1.9 1.9], 0.6, ...
              'Finite slope changes cross NEAR the conjugate point; shrink delta to see them converge.');
P(end+1) = mk('No conjugate point: yp^2 + y^2', 'yp^2 + y^2', 0, 4, 0, 1, [-5 5], 0.5, ...
              'Q > 0: the Jacobi field is sinh(t) and never returns to zero.');
P(end+1) = mk('Arc length (straight lines)', 'sqrt(1+yp^2)', 0, 1, 0, 1, [-5 5], 0.5, ...
              'Straight lines; the Jacobi field is t - a and never returns to zero.');
end

% ---------------------------------------------------------------------------
function s = mk(name, F, a, b, ya, yb, pRange, delta, note)
% MK  One preset record. INPUTS: the fields. OUTPUTS: s struct.
s = struct('name', name, 'F', F, 'a', a, 'b', b, 'ya', ya, 'yb', yb, ...
           'pRange', pRange, 'delta', delta, 'note', note);
end
