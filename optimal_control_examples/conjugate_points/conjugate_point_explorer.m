function app = conjugate_point_explorer(varargin)
%% Purpose:
%
%   An interactive picture of a CONJUGATE POINT. Type a Lagrangian
%   F(t, y, yp), an interval and the boundary values; the explorer finds
%   every extremal (by shooting), and for the one you pick shows five views
%   of the same question -- is this extremal a local minimiser?
%
%     Curves      the extremal and its neighbours: the extremals that leave
%                 the same left point with slope p + delta (a fan at delta,
%                 delta/2, delta/4, delta/8). Where they meet the extremal
%                 again is (in the limit delta -> 0) the conjugate point.
%     Difference  neighbour minus extremal, on the same time axis: its zeros
%                 are the crossings. Tick "divide by delta" and the curves
%                 collapse onto the Jacobi field h (dashed) as delta -> 0.
%     Jacobi      h(t) = dy/dp, the linearised neighbour. Its first zero
%                 after a is the conjugate point t_c.
%     Shooting    r(p) = y(b; p) - yb. Its roots are the extremals; its
%                 slope at a root is h(b).
%     Mode        the lowest eigenfunction eta* of the second variation.
%     Delta J     J[y0 + eps eta*] - J[y0]: a downward parabola means a
%                 direction that LOWERS the cost -- not a minimiser.
%
%   The readout states the Legendre check, t_c against b, the Morse count
%   (negative modes = conjugate points in (a, b)) and the verdict.
%
%   Usage:
%     conjugate_point_explorer                 % opens on the first preset
%     app = conjugate_point_explorer('Visible', 'off');   % scripted/tests
%     app.setPreset(4);  app.selectExtremal(2);  app.setDelta(0.1);
%     disp(app.readout());
%
%% Inputs:
%
%  varargin                 name/value              'Visible' ['on'|'off'],
%                                                   'Preset' [index, 2]
%
%% Outputs:
%
%  app                      struct                  .fig and handles:
%                                                   setPreset(k) solve()
%                                                   selectExtremal(k)
%                                                   setDelta(d) readout()
%                                                   state() shrink()
%                                                   setScaled(tf)
%
%% References:
%
%   [1] I. M. Gelfand and S. V. Fomin, "Calculus of Variations,"
%       Prentice-Hall, 1963, Ch. 5.
%
%% Revision History:
%  M. Casey                                                   (c) 09/24/2026
%  Copyright Coorbital Inc.
%% ------------------------ Begin Code Sequence ---------------------------

ip = inputParser;
ip.addParameter('Visible', 'on');
ip.addParameter('Preset', 2);
ip.parse(varargin{:});
here = fileparts(mfilename('fullpath'));
addpath(here);

presets = cov_presets;
st = struct('prob', [], 'E', [], 'scan', [], 'k', 1, 'S', [], 'V', [], ...
            'fan', [], 'tEnd', [], 'curve', [], 'msg', '', 'diff', []);

% ------------------------------------------------------------ layout
fig = uifigure('Name', 'Conjugate Point Explorer', 'Position', [40 40 1500 920], ...
               'Visible', ip.Results.Visible);
main = uigridlayout(fig, [1 2]);
main.ColumnWidth = {420, '1x'};

ctl = uigridlayout(main, [15 2]);
ctl.RowHeight = {30, 22, 22, 22, 26, 26, 26, 26, 32, 26, 22, 26, 22, 48, '1x'};
ctl.ColumnWidth = {'1x', '1x'};

hdr = uilabel(ctl, 'Text', 'Conjugate Point Explorer', 'FontSize', 18, 'FontWeight', 'bold');
hdr.Layout.Column = [1 2];
uilabel(ctl, 'Text', 'Preset');
ddPreset = uidropdown(ctl, 'Items', [{presets.name}, {'Custom'}], ...
                      'ValueChangedFcn', @(~, ~) onPreset());
ddPreset.Layout.Row = 2;  ddPreset.Layout.Column = 2;
lf = uilabel(ctl, 'Text', 'Lagrangian  F(t, y, yp) =');
lf.Layout.Column = [1 2];
efF = uieditfield(ctl, 'text');
efF.Layout.Column = [1 2];
efA  = numfield(ctl, 'a');    efB  = numfield(ctl, 'b');
efYa = numfield(ctl, 'y(a)'); efYb = numfield(ctl, 'y(b)');
efP1 = numfield(ctl, 'slope scan min'); efP2 = numfield(ctl, 'slope scan max');
btnSolve = uibutton(ctl, 'Text', 'Solve for the extremals', 'FontWeight', 'bold', ...
                    'ButtonPushedFcn', @(~, ~) solve());
btnSolve.Layout.Column = [1 2];
uilabel(ctl, 'Text', 'Extremal');
ddExt = uidropdown(ctl, 'Items', {'(none)'}, 'ValueChangedFcn', @(~, ~) onExtremal());
lblD = uilabel(ctl, 'Text', 'Slope change  delta = 0.5');
lblD.Layout.Column = [1 2];
sld = uislider(ctl, 'Limits', [-1 1], 'Value', 0.5, ...
               'ValueChangedFcn', @(s, ~) setDelta(s.Value));
sld.Layout.Column = [1 2];
cbFan = uicheckbox(ctl, 'Text', 'fan: delta, /2, /4, /8', 'Value', true, ...
                   'ValueChangedFcn', @(~, ~) redrawDelta());
uibutton(ctl, 'Text', 'Shrink delta -> 0', 'ButtonPushedFcn', @(~, ~) shrink());
cbScale = uicheckbox(ctl, 'Text', 'difference plot: divide by delta (-> Jacobi field h)', ...
                     'Value', false, 'ValueChangedFcn', @(~, ~) redrawDelta());
cbScale.Layout.Column = [1 2];
lblNote = uilabel(ctl, 'Text', '', 'WordWrap', 'on', 'FontAngle', 'italic');
lblNote.Layout.Column = [1 2];
txt = uitextarea(ctl, 'Editable', 'off', 'FontName', 'Menlo', 'FontSize', 11, 'WordWrap', 'off');
txt.Layout.Column = [1 2];

right = uigridlayout(main, [4 2]);
right.RowHeight = {'1.15x', '0.8x', '1x', '1x'};
axC = uiaxes(right);  axC.Layout.Column = [1 2];
axD = uiaxes(right);  axD.Layout.Column = [1 2];   % neighbour minus extremal
axH = uiaxes(right);  axR = uiaxes(right);
axM = uiaxes(right);  axJ = uiaxes(right);
for ax = [axC axD axH axR axM axJ]
    ax.Box = 'on';  ax.XGrid = 'on';  ax.YGrid = 'on';  ax.FontSize = 11;
end

app = struct('fig', fig, 'setPreset', @setPreset, 'solve', @solve, ...
             'selectExtremal', @selectExtremal, 'setDelta', @setDelta, ...
             'readout', @() strjoin(txt.Value, newline), 'state', @getState, ...
             'shrink', @shrink, 'setScaled', @setScaled);
setPreset(ip.Results.Preset);

% ------------------------------------------------------------ callbacks
    function s = getState()
        % the LIVE state (an anonymous @() st would freeze it at creation)
        s = st;
    end

    function f = numfield(parent, label)
        % a labelled numeric field occupying one grid cell
        g = uigridlayout(parent, [1 2], 'Padding', 0, 'ColumnSpacing', 4);
        g.ColumnWidth = {'fit', '1x'};
        uilabel(g, 'Text', label);
        f = uieditfield(g, 'numeric');
    end

    function setPreset(k)
        q = presets(k);
        ddPreset.Value = q.name;
        efF.Value = q.F;  efA.Value = q.a;  efB.Value = q.b;
        efYa.Value = q.ya;  efYb.Value = q.yb;
        efP1.Value = q.pRange(1);  efP2.Value = q.pRange(2);
        dmax = max(abs(q.delta)*2, 0.2);
        sld.Limits = [-dmax dmax];  sld.Value = q.delta;
        lblNote.Text = q.note;
        solve();
    end

    function onPreset()
        k = find(strcmp({presets.name}, ddPreset.Value), 1);
        if isempty(k)
            lblNote.Text = 'Custom: edit F, the interval and the boundary values, then Solve.';
        else
            setPreset(k);
        end
    end

    function solve()
        try
            st.prob = cov_problem(efF.Value, efA.Value, efB.Value, efYa.Value, efYb.Value);
            [st.E, st.scan] = cov_extremals(st.prob, [efP1.Value efP2.Value], 300);
        catch err
            st.E = [];  st.msg = err.message;
            txt.Value = {['ERROR: ' err.message]};
            clearAll();
            return
        end
        if isempty(st.E)
            ddExt.Items = {'(none in this slope range)'};
            txt.Value = {sprintf('y'''' = %s', st.prob.gStr), '', ...
                         'No extremal found: r(p) does not change sign in the scanned', ...
                         'slope range. Widen it, or the two-point problem has no solution.'};
            clearAll();
            plotShooting();
            return
        end
        items = arrayfun(@(k) sprintf('%d:  y''(a) = %+.4f,  J = %.5f', k, st.E(k).p, st.E(k).J), ...
                         1:numel(st.E), 'UniformOutput', false);
        ddExt.Items = items;
        selectExtremal(1);
    end

    function onExtremal()
        selectExtremal(find(strcmp(ddExt.Items, ddExt.Value), 1));
    end

    function selectExtremal(k)
        if isempty(st.E), return, end
        st.k = k;
        ddExt.Value = ddExt.Items{k};
        pr = st.prob;
        st.tEnd = pr.b + 0.35*(pr.b - pr.a);
        st.S = cov_shoot(pr, st.E(k).p, st.tEnd);
        st.V = cov_second_variation(pr, st.E(k), 300);
        ys = deval(st.E(k).S.sol, st.V.t);
        epsMax = 0.15*max(1, max(abs(ys(1,:))));
        epsv = linspace(-epsMax, epsMax, 21);
        J0 = cov_functional(pr, st.E(k).S.sol, zeros(size(st.V.eta)));
        dJ = arrayfun(@(e) cov_functional(pr, st.E(k).S.sol, e*st.V.eta) - J0, epsv);
        st.curve = struct('eps', epsv, 'dJ', dJ, ...
                          'quad', st.V.dJquad(end)/st.V.eps(end)^2*epsv.^2);
        plotStatic();
        redrawDelta();
    end

    function setDelta(d)
        if d == 0, d = 1e-6; end
        d = max(min(d, sld.Limits(2)), sld.Limits(1));
        sld.Value = d;
        redrawDelta();
    end

    function setScaled(v)
        % the difference panel: raw (false) or divided by delta (true)
        cbScale.Value = logical(v);
        redrawDelta();
    end

    function shrink()
        d0 = sld.Value;
        for kk = 1:24
            setDelta(d0*0.82^kk);
            drawnow;
            pause(0.04);
        end
    end

% ------------------------------------------------------------ drawing
    function clearAll()
        for ax = [axC axD axH axR axM axJ]
            cla(ax);  legend(ax, 'off');  delete(findall(ax, 'Tag', 'beyondB'));
        end
    end

    function plotStatic()
        pr = st.prob;  S = st.S;  V = st.V;
        tc = S.tConj;
        % Jacobi field
        cla(axH);  hold(axH, 'on');
        tt = linspace(pr.a, S.tStop, 600);
        z = deval(S.sol, tt);
        shadeBeyond(axH, pr.b, S.tStop);
        plot(axH, tt, z(3,:), 'LineWidth', 2, 'Color', [0.10 0.35 0.75]);
        yline(axH, 0, 'k-');
        xline(axH, pr.b, 'k--', 'b', 'LabelVerticalAlignment', 'bottom');
        if ~isempty(tc)
            plot(axH, tc, 0*tc, 'o', 'MarkerSize', 9, 'LineWidth', 2, 'Color', [0 0.6 0]);
        end
        hold(axH, 'off');
        xlim(axH, [pr.a S.tStop]);
        % scale to [a, b]: past b a Jacobi field can grow by orders of magnitude
        hAB = max(abs(z(3, tt <= pr.b)));
        if hAB > 0, ylim(axH, [-1.3 1.3]*hAB); end
        title(axH, 'Jacobi field  h(t) = \partial y/\partial p   (zeros = conjugate points)');
        xlabel(axH, 't');
        % lowest mode
        cla(axM);  hold(axM, 'on');
        plot(axM, V.t, V.eta, 'LineWidth', 2, 'Color', [0.55 0.15 0.65]);
        yline(axM, 0, 'k-');
        tcIn = tc(tc < pr.b);
        for c = tcIn, xline(axM, c, ':', 't_c', 'Color', [0 0.6 0], 'LineWidth', 1.5); end
        hold(axM, 'off');
        xlim(axM, [pr.a pr.b]);
        title(axM, sprintf('Lowest mode \\eta^* of the second variation,  \\lambda_1 = %+.4f', V.lambda(1)));
        xlabel(axM, 't');
        % Delta J along the mode
        cla(axJ);  hold(axJ, 'on');
        if V.lambda(1) < 0, col = [0.80 0.15 0.15]; else, col = [0.10 0.55 0.20]; end
        plot(axJ, st.curve.eps, st.curve.dJ, 'o-', 'LineWidth', 2, 'Color', col, 'MarkerSize', 4);
        plot(axJ, st.curve.eps, st.curve.quad, '--', 'Color', [0.4 0.4 0.4]);
        yline(axJ, 0, 'k-');
        hold(axJ, 'off');
        title(axJ, '\DeltaJ = J[y_0 + \epsilon\eta^*] - J[y_0]   (dashed: \epsilon^2/2 \cdot second variation)');
        xlabel(axJ, '\epsilon');
        plotShooting();
    end

    function plotShooting()
        cla(axR);  hold(axR, 'on');
        plot(axR, st.scan.p, st.scan.r, 'LineWidth', 1.8, 'Color', [0.2 0.2 0.2]);
        yline(axR, 0, 'k-');
        if ~isempty(st.E)
            plot(axR, [st.E.p], 0*[st.E.p], 'o', 'MarkerSize', 8, 'LineWidth', 1.5, 'Color', [0.10 0.35 0.75]);
            e = st.E(st.k);
            plot(axR, e.p, 0, 'o', 'MarkerSize', 11, 'MarkerFaceColor', [0.10 0.35 0.75], 'Color', 'k');
            w = 0.08*diff(st.scan.p([1 end]));
            plot(axR, e.p + [-w w], e.hb*[-w w], '-', 'LineWidth', 2.5, 'Color', [0 0.6 0]);
            text(axR, e.p + w, e.hb*w, sprintf('  slope h(b) = %.3f', e.hb), 'Color', [0 0.5 0]);
        end
        hold(axR, 'off');
        title(axR, 'Shooting function  r(p) = y(b; p) - y_b   (roots = extremals)');
        xlabel(axR, 'initial slope  p = y''(a)');
        fin = st.scan.r(isfinite(st.scan.r));
        if ~isempty(fin)
            lim = prctile(abs(fin), 90);
            if lim > 0, ylim(axR, [-1.2 1.2]*lim); end
        end
    end

    function redrawDelta()
        if isempty(st.E), return, end
        pr = st.prob;  e = st.E(st.k);  d = sld.Value;
        lblD.Text = sprintf('Slope change  delta = %.4g', d);
        if cbFan.Value, ds = d*[1 1/2 1/4 1/8]; else, ds = d; end
        st.fan = arrayfun(@(dd) cov_perturbed(pr, e, dd, st.tEnd), ds);
        % curves
        cla(axC);  hold(axC, 'on');
        shadeBeyond(axC, pr.b, st.tEnd);
        t0 = linspace(pr.a, st.S.tStop, 800);
        z0 = deval(st.S.sol, t0);
        inB = t0 <= pr.b;
        hE = plot(axC, t0(inB), z0(1, inB), 'LineWidth', 3, 'Color', [0.10 0.35 0.75]);
        plot(axC, t0(~inB), z0(1, ~inB), '--', 'LineWidth', 2, 'Color', [0.10 0.35 0.75]);
        cols = [0.95 0.45 0.10; 0.95 0.60 0.25; 0.95 0.72 0.45; 0.95 0.82 0.65];
        hP = gobjects(0);
        for kk = 1:numel(st.fan)
            F1 = st.fan(kk);
            if isempty(F1.S.sol), continue, end
            t1 = linspace(pr.a, F1.S.tStop, 800);
            z1 = deval(F1.S.sol, t1);
            hP(end+1) = plot(axC, t1, z1(1,:), 'LineWidth', 1.6, 'Color', cols(kk,:)); %#ok<AGROW>
            if ~isempty(F1.tCross)
                zc = deval(st.S.sol, F1.tCross(1));
                plot(axC, F1.tCross(1), zc(1), 'o', 'MarkerSize', 9, 'LineWidth', 2, 'Color', [0.8 0.1 0.1]);
            end
        end
        plot(axC, [pr.a pr.b], [pr.ya pr.yb], 's', 'MarkerSize', 9, 'MarkerFaceColor', 'k', 'Color', 'k');
        xline(axC, pr.b, 'k--', 'b', 'LabelVerticalAlignment', 'bottom');
        for c = st.S.tConj
            xline(axC, c, ':', 't_c', 'Color', [0 0.6 0], 'LineWidth', 2);
        end
        hold(axC, 'off');
        xlim(axC, [pr.a st.tEnd]);
        % scale to [a, b] (plus a little beyond): the curves can run away past b
        tView = linspace(pr.a, min(st.tEnd, pr.b + 0.15*(pr.b - pr.a)), 400);
        zv = deval(st.S.sol, tView(tView <= st.S.tStop));
        allY = zv(1,:);
        for kk = 1:numel(st.fan)
            F1 = st.fan(kk);
            if isempty(F1.S.sol), continue, end
            z1 = deval(F1.S.sol, tView(tView <= F1.S.tStop));
            allY = [allY, z1(1,:)]; %#ok<AGROW>
        end
        rngY = [min(allY) max(allY)];
        pad = 0.15*max(diff(rngY), 1e-3);
        ylim(axC, rngY + [-pad pad]);
        if ~isempty(hP)
            legend(axC, [hE hP(1)], {'extremal y_0', 'neighbours (slope p + \delta, ...)'}, ...
                   'Location', 'best');
        end
        title(axC, sprintf('F = %s:  extremal and its neighbours  (red o: first crossing; green: t_c)', pr.Fstr), ...
              'Interpreter', 'none');
        ylabel(axC, 'y');
        drawDifference(cols);
        % mark p + delta on the shooting plot
        plotShooting();
        hold(axR, 'on');
        xline(axR, e.p + d, '-', 'p+\delta', 'Color', [0.95 0.45 0.10], 'LineWidth', 1.5);
        hold(axR, 'off');
        writeReadout();
    end

    function drawDifference(cols)
        % neighbour minus extremal, y(t; p+delta) - y0(t): its zeros ARE the
        % crossings. Divided by delta it tends to the Jacobi field h (dashed)
        % as delta -> 0 -- exactly, for any delta, on a linear problem.
        pr = st.prob;
        scaled = cbScale.Value;
        cla(axD);  hold(axD, 'on');
        shadeBeyond(axD, pr.b, st.tEnd);
        tD = linspace(pr.a, st.tEnd, 800);
        tD = tD(tD <= st.S.tStop);
        z0 = deval(st.S.sol, tD);
        D = nan(numel(st.fan), numel(tD));
        for kk = 1:numel(st.fan)
            F1 = st.fan(kk);
            if isempty(F1.S.sol), continue, end
            ok1 = tD <= F1.S.tStop;
            z1 = deval(F1.S.sol, tD(ok1));
            D(kk, ok1) = z1(1,:) - z0(1, ok1);
            if scaled, D(kk,:) = D(kk,:)/F1.delta; end
            plot(axD, tD, D(kk,:), 'LineWidth', 1.6, 'Color', cols(kk,:));
            if ~isempty(F1.tCross)
                plot(axD, F1.tCross(1), 0, 'o', 'MarkerSize', 9, 'LineWidth', 2, 'Color', [0.8 0.1 0.1]);
            end
        end
        h = z0(3,:);
        if scaled
            plot(axD, tD, h, 'k--', 'LineWidth', 2);
        end
        yline(axD, 0, 'k-');
        xline(axD, pr.b, 'k--', 'b', 'LabelVerticalAlignment', 'bottom');
        for c = st.S.tConj
            xline(axD, c, ':', 't_c', 'Color', [0 0.6 0], 'LineWidth', 2);
        end
        hold(axD, 'off');
        xlim(axD, [pr.a st.tEnd]);
        inView = tD <= min(st.tEnd, pr.b + 0.15*(pr.b - pr.a));
        vals = D(:, inView);
        if scaled, vals = [vals(:); h(inView).']; end
        lim = max(abs(vals(isfinite(vals))));
        if ~isempty(lim) && lim > 0, ylim(axD, [-1.15 1.15]*lim); end
        if scaled
            title(axD, '(y(t; p+\delta) - y_0(t)) / \delta   \rightarrow   Jacobi field h(t)  (dashed)  as \delta \rightarrow 0');
        else
            title(axD, 'Difference  y(t; p+\delta) - y_0(t)   (zeros = crossings)');
        end
        xlabel(axD, 't');
        st.diff = struct('t', tD, 'D', D, 'h', h, 'scaled', scaled);
    end

    function shadeBeyond(ax, b, tEnd)
        % grey band past b; hidden handles survive cla, so delete the old band
        delete(findall(ax, 'Tag', 'beyondB'));
        yl = [-1e6 1e6];
        patch(ax, [b tEnd tEnd b], yl([1 1 2 2]), [0.93 0.93 0.93], 'EdgeColor', 'none', ...
              'HandleVisibility', 'off', 'Tag', 'beyondB');
    end

    function writeReadout()
        pr = st.prob;  e = st.E(st.k);  V = st.V;  tc = st.S.tConj;
        tol = 1e-9*(pr.b - pr.a);
        tcIn = tc(tc < pr.b - tol);  atB = any(abs(tc - pr.b) <= tol);
        legOK = V.Pmin > 0;
        L = {};
        L{end+1} = sprintf('y'''' = %s', pr.gStr);
        L{end+1} = sprintf('extremal %d of %d:  y''(a) = %+.6f', st.k, numel(st.E), e.p);
        L{end+1} = sprintf('                  J     = %.6f', e.J);
        L{end+1} = '';
        L{end+1} = sprintf('LEGENDRE  min F_ypyp = %.4g  %s', V.Pmin, tern(legOK, '> 0  ok', '<= 0  FAILS'));
        if isempty(tc)
            L{end+1} = sprintf('JACOBI    h has no zero in (a, %.3g]', st.S.tStop);
        else
            L{end+1} = sprintf('JACOBI    h = 0 at %s', strjoin(compose('%.6f', tc), ', '));
        end
        L{end+1} = sprintf('          t_c in (a,b): %s', ...
                           tern(isempty(tcIn), 'none', sprintf('%.6f', min(tcIn))));
        if legOK
            L{end+1} = sprintf('MORSE     negative modes %d, conj pts %d  %s', ...
                               V.nNeg, numel(tcIn), tern(V.nNeg == numel(tcIn), 'agree', 'DISAGREE'));
        else
            % the index theorem assumes P > 0; with P < 0 fast wiggles drive
            % the second variation to -inf, so every fine mode is negative
            L{end+1} = sprintf('MORSE     not applicable (needs P > 0); %d negative modes', V.nNeg);
        end
        L{end+1} = sprintf('          lambda_1 = %+.5f', V.lambda(1));
        L{end+1} = sprintf('DELTA J   eps = %.3g on eta*: %+.3e', st.curve.eps(end), st.curve.dJ(end));
        L{end+1} = '';
        L{end+1} = 'NEIGHBOURS  delta     1st crossing  minus t_c';
        for kk = 1:numel(st.fan)
            F1 = st.fan(kk);
            if isempty(F1.tCross)
                L{end+1} = sprintf('            %-9.4g none', F1.delta); %#ok<AGROW>
            elseif isempty(tc)
                L{end+1} = sprintf('            %-9.4g %.6f', F1.delta, F1.tCross(1)); %#ok<AGROW>
            else
                L{end+1} = sprintf('            %-9.4g %-13.6f %+.2e', F1.delta, F1.tCross(1), ...
                                   F1.tCross(1) - tc(1)); %#ok<AGROW>
            end
        end
        L{end+1} = '';
        if V.Pmin < 0
            % Legendre's condition P >= 0 is NECESSARY for a minimum
            L{end+1} = 'VERDICT: F_ypyp < 0 on the extremal: Legendre''s';
            L{end+1} = '         necessary condition fails, NOT a minimiser.';
        elseif ~legOK
            L{end+1} = 'VERDICT: F_ypyp = 0 somewhere: Legendre borderline,';
            L{end+1} = '         no conclusion from these tests.';
        elseif ~isempty(tcIn)
            L{end+1} = sprintf('VERDICT: conjugate point at %.4f < b:', min(tcIn));
            L{end+1} = '         NOT a local minimiser.';
        elseif atB
            L{end+1} = 'VERDICT: conjugate point exactly at b:';
            L{end+1} = '         borderline, no strict conclusion.';
        else
            L{end+1} = 'VERDICT: no conjugate point in (a,b] and';
            L{end+1} = '         Legendre holds: a weak local minimiser.';
        end
        txt.Value = L;
    end
end

% ---------------------------------------------------------------------------
function s = tern(c, a, b)
% TERN  Inline if. INPUTS: c logical, a, b char. OUTPUTS: s char.
if c, s = a; else, s = b; end
end
