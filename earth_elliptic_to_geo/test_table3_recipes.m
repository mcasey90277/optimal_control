% TEST_TABLE3_RECIPES  Unit tests for table3_certified.m and table3_recipes.m
% (pure, no-solve foundation of the Table-3 reproducer engine).
%
% Run:
%   matlab -batch "run('/abs/path/earth_elliptic_to_geo/test_table3_recipes.m')"
here = fileparts(mfilename('fullpath')); cd(here);

c10 = table3_certified(10);
assert(abs(c10.m_f_kg-1377.10)<1e-2 && c10.switches==19 && abs(c10.revs-7.326)<1e-3);
c05 = table3_certified(0.5);
assert(strcmp(c05.anchorSource,'R0law') && abs(c05.tfmin-446.28)<0.1);

r1 = table3_recipes(1);
assert(strcmp(r1.anchor.strategy,'smallN_first') && r1.anchor.nprLo==15 && r1.anchor.nprHi==25);
assert(~isempty(r1.psr) && r1.psr.maxRounds==2);

r05 = table3_recipes(0.5);
assert(strcmp(r05.anchor.strategy,'R0law') && r05.fuel.npr==12 && r05.psr.maxRounds==5);

r10 = table3_recipes(10);
assert(strcmp(r10.anchor.strategy,'coldB') && isempty(r10.psr));

% seeded deep rungs exist but are flagged not-run
r02 = table3_recipes(0.2);  assert(r02.seeded==true);

fprintf('test_table3_recipes PASSED\n');
