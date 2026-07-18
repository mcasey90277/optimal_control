% TEST_VERIFY_ROW  Unit tests for verify_row.m (pure row-vs-certified check).
%
% Run:
%   matlab -batch "run('/abs/path/earth_elliptic_to_geo/test_verify_row.m')"
here = fileparts(mfilename('fullpath')); cd(here);

cert = table3_certified(10);
tol  = struct('m_f_kg',0.5,'revsRel',0.01,'switchesAbs',0);

good = struct('thrustN',10,'m_f_kg',1377.0,'switches',19,'revs',7.33);
assert(verify_row(good,cert,tol)==true);

bad  = struct('thrustN',10,'m_f_kg',1300,'switches',19,'revs',7.33);
threw=false; try, verify_row(bad,cert,tol); catch, threw=true; end
assert(threw,'verify_row must throw on a mass mismatch');

fprintf('test_verify_row PASSED\n');
