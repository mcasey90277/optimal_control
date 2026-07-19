function recipe = table3_recipes(thrustN)
% TABLE3_RECIPES  Per-rung recipe registry for the Table-3 reproducer engine.
%
% Pure lookup: no file I/O, no solves, no side effects. Each rung's recipe
% captures the exact proven knobs (harvested from the campaign's hand
% scripts, see process/CAMPAIGN.md and process/DESIGN_thrust_ladder.md) needed by
% reproduce_row.m to re-solve that row from scratch: which anchor strategy
% to use (coldB/chain/smallN_first/R0law), the fuel-stage node count and
% warm-start source, and an optional PSR (post-solution refinement) pass.
% 0.2 N and 0.1 N are seeded (recipe present, `.seeded==true`) but have not
% been run to a certified row in this build.
%
% Every anchor struct carries the full field set
% (strategy/npr/nprLo/nprHi/mtMaxIter/warmFrom) with unused fields set to
% [], so callers never need an isfield guard. Likewise every non-empty psr
% struct carries the full field set
% (maxRounds/nbr/globalEvery/globalFactor).
%
% INPUTS:
%   thrustN - max thrust level [N]; one of 10, 5, 2.5, 1, 0.5, 0.2, 0.1    [scalar]
%
% OUTPUTS:
%   recipe - struct with fields:
%     .thrustN     - echoes the input thrust level [N]                        [scalar]
%     .anchor      - struct: .strategy ('coldB'|'chain'|'smallN_first'|
%                    'R0law'), .npr, .nprLo, .nprHi, .mtMaxIter, .warmFrom
%                    (previous-rung thrust [N] to warm-start from, or []
%                    for a cold/no-chain anchor)                               [struct]
%     .fuel        - struct: .npr, .seedThr, .maxIter, .warmFrom (previous-
%                    rung thrust [N] to warm-start the fuel stage from, or
%                    [] for none)                                              [struct]
%     .psr         - struct: .maxRounds, .nbr, .globalEvery, .globalFactor
%                    (post-solution refinement pass), or [] if this rung
%                    needs no PSR                                          [struct|[]]
%     .tfmin_or_R0 - R0-law tfmin estimate R0const/thrustN [ND] for rungs
%                    whose anchor strategy is 'R0law'; [] otherwise           [scalar]
%     .seeded      - true if this recipe is a not-yet-executed deep-rung
%                    seed (0.2/0.1 N); false for the five proven/attained
%                    rungs (10/5/2.5/1/0.5 N)                                [logical]
%
% REFERENCES:
%   [1] Haberkorn, Martinon, Gergaud, "Low Thrust Minimum-Fuel Orbital
%       Transfer: A Homotopic Approach," JGCD 27(6), 2004, Table 3.
%   [2] docs/superpowers/specs/2026-07-18-table3-reproducer-engine-design.md
%       Sec. 4 (anchor strategies) and Sec. 6 (recipe registry table).
%   [3] earth_elliptic_to_geo/process/CAMPAIGN.md (harvested per-rung numbers).

if nargin < 1
    error('table3_recipes:badInput', 'thrustN is required');
end

R0const = 223.14;  % single R0-law constant behind every R0law anchor estimate

tol = 1e-9;
rungs = { ...
    10,  make_recipe(10,  make_anchor('coldB',        25, [],  [],  [],  []  ), ...
                          make_fuel(25, 0.4, 1500, []), ...
                          [], false); ...
    5,   make_recipe(5,   make_anchor('chain',        25, [],  [],  300, 10  ), ...
                          make_fuel(25, 0.4, 1500, 10), ...
                          [], false); ...
    2.5, make_recipe(2.5, make_anchor('chain',        25, [],  [],  300, 5   ), ...
                          make_fuel(25, 0.4, 1500, 5), ...
                          [], false); ...
    1,   make_recipe(1,   make_anchor('smallN_first', [], 15,  25,  [],  2.5 ), ...
                          make_fuel(25, 0.4, 1500, 2.5), ...
                          make_psr(2, 2, [], []), false); ...
    0.5, make_recipe(0.5, make_anchor('R0law',         [], [],  [],  [],  1   ), ...
                          make_fuel(12, 0.4, 1500, 1), ...
                          make_psr(5, 2, 3, 1.3), false); ...
    0.2, make_recipe(0.2, make_anchor('chain',         12, [],  [],  [],  0.5 ), ...
                          make_fuel(10, 0.4, 1500, 0.5), ...
                          make_psr(5, 2, 3, 1.3), true); ...
    0.1, make_recipe(0.1, make_anchor('chain',         12, [],  [],  [],  0.2 ), ...
                          make_fuel(8, 0.4, 1500, 0.2), ...
                          make_psr(6, 2, 3, 1.3), true); ...
};

for k = 1:size(rungs,1)
    if abs(rungs{k,1} - thrustN) < tol
        recipe = rungs{k,2};
        if strcmp(recipe.anchor.strategy, 'R0law')
            recipe.tfmin_or_R0 = R0const / thrustN;
        else
            recipe.tfmin_or_R0 = [];
        end
        return;
    end
end

error('table3_recipes:unknownThrust', ...
    'No recipe registered for thrustN=%g (known: 10, 5, 2.5, 1, 0.5, 0.2, 0.1)', thrustN);

end

function anchor = make_anchor(strategy, npr, nprLo, nprHi, mtMaxIter, warmFrom)
% MAKE_ANCHOR  Build one anchor-stage recipe struct with the full,
% isfield-guard-free field set (unused fields set to []).
%
% INPUTS:
%   strategy  - 'coldB'|'chain'|'smallN_first'|'R0law'                    [char]
%   npr       - nodes-per-rev for a single-grid anchor, or []             [scalar|[]]
%   nprLo     - low-grid nodes-per-rev (smallN_first stage 1), or []      [scalar|[]]
%   nprHi     - high-grid nodes-per-rev (smallN_first stage 2), or []     [scalar|[]]
%   mtMaxIter - min-time solver max iterations for a 'chain' anchor, or []  [scalar|[]]
%   warmFrom  - previous-rung thrust [N] to warm-start from, or []        [scalar|[]]
%
% OUTPUTS:
%   anchor - struct with fields .strategy/.npr/.nprLo/.nprHi/.mtMaxIter/.warmFrom
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-18-table3-reproducer-engine-design.md Sec. 4.
anchor = struct('strategy', strategy, 'npr', npr, 'nprLo', nprLo, ...
    'nprHi', nprHi, 'mtMaxIter', mtMaxIter, 'warmFrom', warmFrom);
end

function fuel = make_fuel(npr, seedThr, maxIter, warmFrom)
% MAKE_FUEL  Build one fuel-stage recipe struct.
%
% INPUTS:
%   npr      - nodes-per-rev for the fuel-stage collocation grid          [scalar]
%   seedThr  - homotopy seed threshold                                    [scalar]
%   maxIter  - fuel-stage solver max iterations                           [scalar]
%   warmFrom - previous-rung thrust [N] to warm-start the fuel stage
%              from, or [] for none                                       [scalar|[]]
%
% OUTPUTS:
%   fuel - struct with fields .npr/.seedThr/.maxIter/.warmFrom
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-18-table3-reproducer-engine-design.md Sec. 6.
fuel = struct('npr', npr, 'seedThr', seedThr, 'maxIter', maxIter, 'warmFrom', warmFrom);
end

function psr = make_psr(maxRounds, nbr, globalEvery, globalFactor)
% MAKE_PSR  Build one post-solution-refinement (PSR) recipe struct.
%
% INPUTS:
%   maxRounds    - maximum PSR refinement rounds                          [scalar]
%   nbr          - neighbor-basin count per round                         [scalar]
%   globalEvery  - run a global (vs. local) round every N rounds, or []   [scalar|[]]
%   globalFactor - global-round step-size multiplier, or []               [scalar|[]]
%
% OUTPUTS:
%   psr - struct with fields .maxRounds/.nbr/.globalEvery/.globalFactor
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-18-table3-reproducer-engine-design.md Sec. 6.
psr = struct('maxRounds', maxRounds, 'nbr', nbr, ...
    'globalEvery', globalEvery, 'globalFactor', globalFactor);
end

function recipe = make_recipe(thrustN, anchor, fuel, psr, seeded)
% MAKE_RECIPE  Assemble one rung's full recipe struct.
%
% INPUTS:
%   thrustN - max thrust level [N]                                        [scalar]
%   anchor  - anchor-stage struct (see make_anchor)                       [struct]
%   fuel    - fuel-stage struct (see make_fuel)                           [struct]
%   psr     - PSR struct (see make_psr), or [] if this rung needs no PSR  [struct|[]]
%   seeded  - true if this is a not-yet-executed deep-rung seed           [logical]
%
% OUTPUTS:
%   recipe - struct with fields .thrustN/.anchor/.fuel/.psr/.seeded
%            (.tfmin_or_R0 is filled in by the caller, table3_recipes)
%
% REFERENCES:
%   [1] docs/superpowers/specs/2026-07-18-table3-reproducer-engine-design.md Sec. 6.
recipe = struct('thrustN', thrustN, 'anchor', anchor, 'fuel', fuel, ...
    'psr', psr, 'seeded', seeded);
end
