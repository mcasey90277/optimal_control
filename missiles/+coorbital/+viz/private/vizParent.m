function [hFig,hAx] = vizParent(parent,nAx,figName,visName)
%% Purpose:
%
%  Resolve the 'Parent' option shared by every function in coorbital.viz into
%  a figure handle and the axes to draw into. This is the single place that
%  decides whether a new figure is created, so no plotting routine in the
%  package can create one behind a caller's back.
%
%  Three forms of parent are accepted:
%
%    []                 Create a new figure and nAx axes inside it. This is
%                       the only branch that calls figure().
%
%    axes array         Draw into exactly these axes. numel must equal nAx.
%                       Nothing is created, nothing is cleared, and the
%                       axes' own visibility is left alone.
%
%    figure / panel     Create nAx axes inside the given container. The
%                       container itself is not created and its visibility
%                       is left alone.
%
%  Axes are created with an explicit OuterPosition rather than through
%  subplot, because subplot DELETES any existing axes it overlaps. That is
%  the opposite of composing into a caller's figure, which is the whole point
%  of the option.
%
%% Inputs:
%
%  parent           [] / graphics handle        Empty for a new figure; an
%                                               [1 x nAx] axes array to draw
%                                               into; or a scalar figure or
%                                               panel to create axes inside
%
%  nAx              [1 x 1]                     Number of axes the caller
%                                               needs (-)
%
%  figName          Char [1 x n]                Figure Name property, used
%                                               only when a figure is created
%
%  visName          Char [1 x m]                'on' or 'off', the created
%                                               figure's Visible property.
%                                               Used only when a figure is
%                                               created; a supplied parent
%                                               keeps whatever it has
%
%% Outputs:
%
%  hFig             [1 x 1] figure              The figure the axes live in,
%                                               created or inherited
%
%  hAx              [1 x nAx] axes              Axes to draw into, ordered top
%                                               to bottom when created here
%
%% Notes:
%
%  A private function, so it carries no self-demo: nothing outside
%  +coorbital/+viz can call it, and inside the package it is reached by plain
%  name rather than by namespace.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

    assert(nAx >= 1 && nAx == fix(nAx),'nAx must be a positive whole number.');
    assert(any(strcmp(visName,{'on','off'})), ...
        'Visible must be ''on'' or ''off'', not ''%s''.',visName);

%% No parent: the one branch allowed to create a figure:
    if isempty(parent)
              hFig = figure('color',[1 1 1],'name',figName, ...
                            'numbertitle','off','visible',visName);
               hAx = makeAxes(hFig,nAx);
        return;
    end

    assert(all(isgraphics(parent(:))), ...
        'Parent must be a valid graphics handle, an axes array, or [].');

%% An axes array: compose into it, create nothing:
    if all(isgraphics(parent(:),'axes'))
        assert(numel(parent) == nAx, ...
            ['Parent carries %d axes; this plot needs exactly %d. Pass that ' ...
             'many axes, or pass the figure and let the axes be created.'], ...
            numel(parent),nAx);
               hAx = reshape(parent,1,nAx);
              hFig = ancestor(hAx(1),'figure');
        assert(~isempty(hFig),'the supplied axes are not parented to a figure.');
        return;
    end

%% A container: create the axes inside it, but not the container:
    assert(isscalar(parent), ...
        'Parent must be a single container or an array of axes, not a mix.');
              hFig = ancestor(parent,'figure');
    assert(~isempty(hFig), ...
        'Parent is not a figure and is not parented to one.');
               hAx = makeAxes(parent,nAx);
end

function hAx = makeAxes(container,nAx)
%% Purpose:
%
%  Create nAx axes stacked top to bottom inside a container, each occupying
%  an equal share of its height.
%
%% Inputs:
%
%  container        [1 x 1] graphics            Figure or panel to parent the
%                                               axes to
%
%  nAx              [1 x 1]                     Number of axes to create (-)
%
%% Outputs:
%
%  hAx              [1 x nAx] axes              The created axes, first at the
%                                               top
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

               hAx = gobjects(1,nAx);
    for ka = 1:nAx
              yLow = (nAx - ka)./nAx;
           hAx(ka) = axes('Parent',container, ...
                          'OuterPosition',[0 yLow 1 1./nAx]);
    end
end
