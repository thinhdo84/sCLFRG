% Calculate jacobian matrix
function [A, B] = linearize_cstr(xe,ue)
q =100;
CA0 = 1;
T0 = 350; Tc0 = 350;
V = 100;
hA = 7e5;
k0 = 7.2e10;
EsR =1e4; %E/R
DeltaH = 2e5;
rho = 1e3; rhoc = rho;
Cp = 1; Cpc = 1;

syms CA T qc
f = [q*(CA0-CA)/V - k0*CA*exp(-EsR/T);
 q*(T0 - T)/V - (-DeltaH)*k0*CA*exp(-EsR/T)/(rho*Cp)+...
    (rhoc*Cpc/(rho*Cp*V))*qc*(1-exp(-hA/(qc*rho*Cp)))*(Tc0-T) ;];

A=jacobian(f,[CA; T]);
A = eval(subs(A, {CA,T, qc}, {xe(1), xe(2), ue}));

B=jacobian(f,qc);
B = eval(subs(B, {CA,T, qc}, {xe(1), xe(2), ue}));



end