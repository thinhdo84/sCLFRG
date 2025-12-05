function xd = bioreact(x,u, para)
% Created by Huu-Thinh DO, 2024
% Based on:"Multistep Nonlinear Predictive Controller
% David D. Brengel and Warren D. Seider, 1989.
% eq (28)-(29)
% This vector field describes the deviation dynamics around a working
% equilibrium point of a bioreactor process. The objective is to drive
%(x,u) to the origin.

% x : the system's state in R2,
% u : the control in R,
% para: the structure containing all constant and operating point.
% output xd: the first derivative of the state x.

C1S=para.C1S;C2S=para.C2S;DaS=para.DaS;
beta=para.beta;
gamma=para.gamma;

xd =[-(C1S + x(1))*(1 + u)+...
DaS*(1 - x(2)- C2S)*exp((C2S + x(2))/gamma)*(C1S+x(1));
-(C2S + x(2))*(1 + u)+...
DaS*(1 - x(2)- C2S)*exp((C2S + x(2))/gamma)*(1+beta)*(C1S+x(1))/(1+beta-x(2)-C2S)];
end