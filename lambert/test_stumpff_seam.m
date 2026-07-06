zL = 1e-4 - 1e-12;
zR = 1e-4 + 1e-12;
[CL, SL] = stumpff(zL);
[CR, SR] = stumpff(zR);
disp([abs(CL-CR), abs(SL-SR)]);
