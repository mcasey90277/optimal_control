% test_gergaud_row.m
here = fileparts(mfilename('fullpath')); cd(here);
inp = struct('thrustN',10,'tfmin_ND',22.2206,'ctf',1.5,'tf_ND',33.331, ...
    'm_f_kg',1377.10,'switches',19,'revs',7.326,'edge',0.999,'incl_deg',0.0, ...
    'defect',6e-15,'certified',true,'note','');
row = gergaud_row(inp);
assert(abs(row.prop_kg-(1500-1377.10)) < 1e-6, 'prop = m0 - m_f');
assert(row.revs_paper==7.5, 'paper revs lookup 10 N -> 7.5');
assert(row.dV_kms > 0, 'dV positive');
s = gergaud_row_str(row);
assert(ischar(s) && contains(s,'1377.10'), 'row string carries m_f');
% uncertified row must be flagged, not silently formatted as certified
inp2 = inp; inp2.thrustN=0.2; inp2.certified=false; inp2.note='not attained (0.5 N wall)';
s2 = gergaud_row_str(gergaud_row(inp2));
assert(contains(lower(s2),'not attained') || contains(s2,'UNCERTIFIED'), 'uncertified flagged');
% FIX I-2: a CERTIFIED row with a non-empty note (e.g. the 0.5 N anchor-
% free tfmin estimate, or the I-1 PSR-skipped-for-custom-endpoints note)
% must still print that note -- a printed row must never be mistaken for
% a fully-clean certified one just because certified==true.
inp3 = inp; inp3.certified=true; ...
    inp3.note='PSR switch-refinement skipped for custom endpoints (research-probe): reported solution is the un-refined fuel solve';
s3 = gergaud_row_str(gergaud_row(inp3));
assert(contains(s3, inp3.note), 'certified-but-footnoted row must still print its note');
assert(~contains(s3, 'UNCERTIFIED'), 'certified row must not carry the UNCERTIFIED banner');
fprintf('test_gergaud_row PASSED\n');
