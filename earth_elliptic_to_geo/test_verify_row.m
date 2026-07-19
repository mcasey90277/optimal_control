% TEST_VERIFY_ROW  One-sided verify_row: reproduced mass must be at least the
% campaign floor; higher mass (better) always passes and is flagged improved;
% structure (switches/revs) is reported, not gated.
here = fileparts(mfilename('fullpath')); cd(here);
cert = table3_certified(10);          % 1377.10 kg / 19 sw / 7.326 rev
tol  = struct('m_f_kg', 0.5);

% (a) exactly at the campaign mass -> passes, not flagged improved
atCampaign = struct('thrustN',10,'m_f_kg',cert.m_f_kg,'switches',19,'revs',7.326);
[p, info] = verify_row(atCampaign, cert, tol);
assert(p == true && info.improved == false, 'at-campaign row should pass, not improved');

% (b) BETTER mass with a DIFFERENT structure (18 sw / 7.56 rev) -> passes,
%     flagged improved, structure NOT gated
better = struct('thrustN',10,'m_f_kg',1378.46,'switches',18,'revs',7.560);
[p2, info2] = verify_row(better, cert, tol);
assert(p2 == true, 'a higher-mass row must pass regardless of structure');
assert(info2.improved == true, 'a higher-mass row must be flagged improved');
assert(abs(info2.improvedKg - (1378.46 - cert.m_f_kg)) < 1e-6, 'improvedKg correct');

% (c) mass BELOW the floor -> throws worseThanCampaign
worse = struct('thrustN',10,'m_f_kg',cert.m_f_kg - 1.0,'switches',19,'revs',7.326);
threw = false;
try
    verify_row(worse, cert, tol);
catch ME
    threw = strcmp(ME.identifier, 'verify_row:worseThanCampaign');
end
assert(threw, 'a mass below the floor must throw verify_row:worseThanCampaign');

fprintf('test_verify_row PASSED\n');
