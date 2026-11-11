function st = pole_series_stats(tDays, gdop, tLastDays, dtSec)
%% Purpose:
%
%   SUMMARY STATISTICS of the study's south-pole GDOP series (no verdicts).
%   The service window is every sample with tDays >= tLastDays. A sample
%   counts as SERVED iff GDOP < 6; NaN (fewer than four in view), Inf and
%   GDOP >= 6 (including exactly 6) are outages. Median and 95th percentile
%   are over the FINITE service samples (NaN when there are none).
%
%     firstFiniteDay  first tDays with finite GDOP (whole series; NaN if never)
%     firstLt6Day     first tDays with GDOP < 6 (whole series; NaN if never)
%     svcSamples      service-window samples
%     svcDays         svcSamples*dtSec/86400
%     fracLt6         fraction of service samples with GDOP < 6
%     fracFinite      fraction of service samples with finite GDOP
%     medianGdop      median of the finite service samples
%     p95Gdop         95th percentile of the finite service samples
%     longestOutageH  longest run of consecutive service-window outages
%                     times dtSec, in hours
%     n6to7           service samples with 6 <= GDOP < 7
%     nLt6            service samples with GDOP < 6
%
%   Self-demo (nargin == 0): a hand-built series with one 3-sample outage.
%
%% Inputs:
%
%  tDays                    [N x 1]                 Epochs [d]
%
%  gdop                     [N x 1]                 Pole GDOP (NaN, Inf allowed)
%
%  tLastDays                scalar                  Service-window start [d]
%
%  dtSec                    scalar                  Sample spacing [s]
%
%% Outputs:
%
%  st                       struct                  The fields listed above
%                                                   (scalars)
%
%% Revision History:
%  M. Casey                                                   (c) 11/10/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

if nargin == 0
    tD = (0:9).'*300/86400;
    g  = [NaN NaN 5 5 7 NaN 6 5 4 5].';
    st = pole_series_stats(tD, g, tD(3), 300);
    fprintf(['pole_series_stats demo: first finite %.4f d, service %d samples, ' ...
             'GDOP < 6 %.3f, longest outage %.2f h, 6 <= GDOP < 7: %d\n'], ...
            st.firstFiniteDay, st.svcSamples, st.fracLt6, st.longestOutageH, st.n6to7);
    return
end
if numel(tDays) ~= numel(gdop)
    error('pole_series_stats:shape', 'tDays (%d) and gdop (%d) differ in length', ...
          numel(tDays), numel(gdop));
end
tDays = tDays(:);  gd = gdop(:);
svc = tDays >= tLastDays;
gs  = gd(svc);
ok6 = gs < 6;                                          % NaN and Inf count as outages
gf  = gs(isfinite(gs));
st.firstFiniteDay = first_day(tDays, isfinite(gd));
st.firstLt6Day    = first_day(tDays, gd < 6);
st.svcSamples     = numel(gs);
st.svcDays        = numel(gs)*dtSec/86400;
st.fracLt6        = mean(ok6);
st.fracFinite     = mean(isfinite(gs));
if isempty(gf)
    st.medianGdop = NaN;  st.p95Gdop = NaN;
else
    st.medianGdop = median(gf);
    st.p95Gdop    = prctile(gf, 95);
end
st.longestOutageH = longest_run(~ok6)*dtSec/3600;
st.n6to7          = nnz(gs >= 6 & gs < 7);
st.nLt6           = nnz(ok6);
end

%% ======================================================================
function d = first_day(tDays, mask)
% FIRST_DAY  First time where mask holds, NaN if never.
%   INPUTS:  tDays [N x 1]; mask [N x 1] logical.
%   OUTPUTS: d scalar [d].
k = find(mask, 1);
if isempty(k), d = NaN; else, d = tDays(k); end
end

function n = longest_run(mask)
% LONGEST_RUN  Length of the longest run of true samples.
%   INPUTS:  mask [N x 1] logical.
%   OUTPUTS: n scalar (samples).
dm = diff([false; mask(:); false]);
n  = max([0; find(dm == -1) - find(dm == 1)]);
end
