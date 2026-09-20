function [lo, k] = between_sample_bound(t, v, L)
%% Purpose:
%
%   A LOWER-BOUND ESTIMATE of min v(t) over the whole interval [t(1), t(end)], from
%   samples of v and a bound L on its slope. A sampled minimum only says what
%   v was AT the samples; a hypothesis of the form "v(t) > 0 for all t" needs
%   the gaps between them covered too.
%
%   On one interval [t_k, t_k+1] with |v'| <= L, the function can fall from
%   v_k at rate L and must climb back to v_k+1 at rate L, so
%
%       min v  >=  (v_k + v_k+1)/2  -  L (t_k+1 - t_k)/2 ,
%
%   and the bound returned is the smallest of these over all intervals.
%
%  ASSUMPTIONS / NOTES:
%
% • L(k) is a bound on |v'| AT sample k; per interval the LARGER of its two
%   ends is used. That covers the interval exactly when the true slope bound
%   does not exceed its end values inside it -- an ESTIMATE on a finely
%   sampled flight, not a validated enclosure. It is never looser than
%   L dt/2 below the sampled minimum, and it closes on the true minimum as
%   the sampling is refined (tests/test_between_sample_bound).
% • Use: the minimum-time hypotheses H2 (|lam_v| > 0) and H3 (Q_mt > 0) in
%   mintime_hypothesis_gates, where one quantity bounds both slopes:
%   |d|lam_v|/dt| <= |lam_r| in the CR3BP (Coriolis is skew), and dQ_mt/dt =
%   (d|lam_v|/dt)/m exactly on an all-burn arc.
%
%% Inputs:
%
%  t                        [n x 1]                 sample times, increasing
%  v                        [n x 1]                 the samples
%  L                        [n x 1]                 slope bound at each sample
%
%% Outputs:
%
%  lo                       double                  lower bound on min v
%  k                        int                     the interval [t(k), t(k+1)]
%                                                   that gives it
%
%% Revision History:
%  M. Casey                                                   (c) 09/18/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

t = t(:);  v = v(:);  L = L(:);
assert(numel(t) >= 2 && numel(v) == numel(t) && numel(L) == numel(t), 'between_sample_bound:size', ...
       't, v and L must be vectors of one length (at least two samples)');
assert(isreal(t) && isreal(v) && isreal(L) && all(isfinite(t)) && all(isfinite(v)) && all(isfinite(L)), ...
       'between_sample_bound:finite', 'times, samples and slope bounds must be real and finite');
assert(all(L >= 0), 'between_sample_bound:slope', 'a slope bound is a magnitude: it cannot be negative');
dt = diff(t);
assert(all(dt > 0), 'between_sample_bound:time', 'the sample times must be strictly increasing');
perInterval = 0.5*(v(1:end-1) + v(2:end)) - 0.5*max(L(1:end-1), L(2:end)).*dt;
[lo, k] = min(perInterval);
end
