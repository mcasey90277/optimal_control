% TEST_FOC_MANIFEST  Field-contract test for all campaign manifests.
root = fileparts(fileparts(mfilename('fullpath'))); cd(root); setup_verify_common;
names = {'earth_mee','earth_cr3bp','tulip','elfo_fuel','elfo_mintime','toy'};
need  = {'name','nx','nu','dirRows','thrRow','massRow','timeRow', ...
         'autonomous','horizonKind','massFreeAtTf'};
for q = 1:numel(names)
    m = foc_manifest(names{q});
    for f = 1:numel(need), assert(isfield(m, need{f}), '%s missing %s', names{q}, need{f}); end
end
m = foc_manifest('earth_mee');   assert(m.nx==7 && m.massRow==6 && m.timeRow==7 && m.autonomous);
m = foc_manifest('earth_cr3bp'); assert(m.nx==7 && ~m.autonomous);
m = foc_manifest('tulip');       assert(m.nx==8 && m.massRow==7 && m.timeRow==8);
m = foc_manifest('elfo_fuel');   assert(m.nx==9 && strcmp(m.horizonKind,'freetf-cscale'));
m = foc_manifest('elfo_mintime');assert(m.nu==3 && isempty(m.thrRow));
fprintf('test_foc_manifest: ALL PASS\n');
