% TEST_CR3BP_FP  Unit tests: fingerprint build/check + thrust_tag.
here = fileparts(mfilename('fullpath'));  addpath(here);
p  = cr3bp_lt_params(0.025, 15, 2100);
fp = cr3bp_fingerprint(p, struct('tf', 7.23, 'insertion', 'campaign'));
assert(fp.thrustN==0.025 && fp.m0kg==15 && fp.ispS==2100, 'core fields');
assert(abs(fp.Tmax - p.Tmax) < 1e-15 && abs(fp.muStar - p.muStar) < 1e-15, 'derived fields');
assert(fp.tf==7.23 && strcmp(fp.insertion,'campaign'), 'extra fields merged');
% (a) match -> silent
S1 = struct('fp', fp);
check_cr3bp_fp(S1, fp, 'file.mat', 'tag');
% (b) legacy (no fp) -> warn, not error
wOld = warning('off','all'); lastwarn('');
check_cr3bp_fp(struct('X',1), fp, 'file.mat', 'tag');
[~, wid] = lastwarn; assert(strcmp(wid,'check_cr3bp_fp:noFingerprint'), 'legacy warn');
% (c) schema-older (cached fp missing a new field) -> warn
fpOld = rmfield(fp, 'insertion');  lastwarn('');
check_cr3bp_fp(struct('fp',fpOld), fp, 'file.mat', 'tag');
[~, wid] = lastwarn; assert(strcmp(wid,'check_cr3bp_fp:schemaOlder'), 'schema warn');
warning(wOld);
% (d) mismatch -> hard error naming the field
fpBad = fp;  fpBad.thrustN = 0.020;
ok = false;
try, check_cr3bp_fp(struct('fp',fpBad), fp, 'file.mat', 'tag');
catch err, ok = strcmp(err.identifier,'check_cr3bp_fp:mismatch') && contains(err.message,'thrustN'); end
assert(ok, 'mismatch must hard-error naming the field');
% (e) thrust_tag
assert(isempty(thrust_tag(0.025)), 'nominal tag must be empty');
assert(strcmp(thrust_tag(0.020), '_T20mN'), '20 mN tag');
assert(strcmp(thrust_tag(0.0325), '_T32p5mN'), 'fractional-mN tag');
fprintf('test_cr3bp_fp: ALL PASS\n');
