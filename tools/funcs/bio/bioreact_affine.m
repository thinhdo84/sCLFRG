function [f, g] = bioreact_affine(x,para)
C1S=para.C1S;C2S=para.C2S;DaS=para.DaS;
beta=para.beta;
gamma=para.gamma;

f = [     - C1S - x(1) - DaS*exp((C2S + x(2))/gamma)*(C1S + x(1))*(C2S + x(2) - 1);
    (DaS*exp((C2S + x(2))/gamma)*(C1S + x(1))*(beta + 1)*(C2S + x(2) - 1))/(C2S - beta + x(2) - 1) - x(2) - C2S];
g = [- C1S - x(1);
    - C2S - x(2)];