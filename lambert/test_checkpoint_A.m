[C0,S0] = stumpff(0)
[Cp,Sp] = stumpff(pi^2);  [Cp - 2/pi^2,  Sp - 1/pi^2]
[Cn,Sn] = stumpff(-4);
[Cn - (cosh(2)-1)/4,  Sn - (sinh(2)-2)/8]
zL = 1e-4-1e-12; zR = 1e-4+1e-12;
[CL,SL] = stumpff(zL); [CR,SR] = stumpff(zR);
[abs(CL-CR), abs(SL-SR)]
