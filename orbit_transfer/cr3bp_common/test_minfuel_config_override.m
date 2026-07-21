% TEST_MINFUEL_CONFIG_OVERRIDE  Override arg + nominal invariance.
here = fileparts(mfilename('fullpath'));  addpath(here);
c0 = minfuel_config();
c1 = minfuel_config(struct('thrustN', 0.020));
assert(c1.thrustN == 0.020, 'override applied');
c1b = rmfield(c1, 'thrustN');  c0b = rmfield(c0, 'thrustN');
assert(isequal(strip_handles(c1b), strip_handles(c0b)), 'override must change ONLY the named field');
assert(strcmp(func2str(c1.fname), func2str(c0.fname)) && strcmp(func2str(c1.fparse), func2str(c0.fparse)), 'handle sources unchanged');
c2 = minfuel_config(struct());
assert(isequal(strip_handles(c2), strip_handles(c0)), 'empty override = default');
assert(isequal(strip_handles(minfuel_config()), strip_handles(c0)), 'no-arg call unchanged');
fprintf('test_minfuel_config_override: ALL PASS\n');

function s = strip_handles(s)
fn = fieldnames(s);
for k = 1:numel(fn)
    if isa(s.(fn{k}), 'function_handle'), s = rmfield(s, fn{k}); end
end
end
