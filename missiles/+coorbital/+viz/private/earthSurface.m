function hSurf = earthSurface(hAx,nCell)
%% Purpose:
%
%  Draw the reference sphere this library flies over into an existing axes,
%  in kilometres, and hand back the surface handle. Factored out of
%  coorbital.viz.globe3D so that the animated globe can render the SAME body
%  rather than growing a second copy of it that then drifts out of step. It
%  draws the planet and nothing else: no view, no lighting, no limits, no
%  labels. Those belong to the caller, because a still figure and a rotating
%  movie want different ones.
%
%  The radius is coorbital.util.missileConst().rE, converted to kilometres
%  once, here. Nothing in the package may draw a sphere of its own radius.
%
%% Inputs:
%
%  hAx              [1 x 1] axes                Axes to draw into. Its hold
%                                               state is left as found
%
%  nCell            [1 x 1]                     Mesh resolution: the sphere is
%                                               built on an nCell x nCell grid
%                                               of quadrilaterals (-). The mesh
%                                               edges are the graticule, so
%                                               this also sets its spacing:
%                                               nCell = 24 gives 15 deg cells
%
%% Outputs:
%
%  hSurf            [1 x 1] surface             The drawn sphere, radius rE in
%                                               km, Tag 'earthSurface'
%
%% Notes:
%
%  Opaque faces and no light source. A lit sphere reads better on screen but
%  needs a renderer that headless MATLAB cannot be relied on to have; the
%  graticule drawn by the mesh edges carries the three-dimensional shape
%  instead, and costs nothing.
%
%  Built with the low-level surface() rather than surf(). surf() runs newplot,
%  which CLEARS the axes it is handed; that would silently wipe whatever a
%  caller composing into an existing axes had already drawn there. surface()
%  adds a child and touches nothing else.
%
%  A private function, so it carries no self-demo: nothing outside
%  +coorbital/+viz can call it.
%
%% Revision History:
%  Michael Casey                                                08/07/2026
%  Copyright 2026 Coorbital, Inc.
%% ------------------------ Begin Code Sequence ---------------------------

                 c = coorbital.util.missileConst();
             rEkm  = c.rE./1000;

%% Unit sphere on the requested mesh, scaled to the planet radius in km:
        [xU,yU,zU] = sphere(nCell);
             hSurf = surface(hAx,rEkm.*xU,rEkm.*yU,rEkm.*zU, ...
                             'FaceColor',[0.84 0.88 0.93], ...
                             'EdgeColor',[0.55 0.63 0.73], ...
                             'EdgeAlpha',0.45, ...
                             'FaceLighting','none', ...
                             'Tag','earthSurface');
end
