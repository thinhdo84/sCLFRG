function [A, B] = linearize_Bioreact(para)
% Created by Huu-Thinh DO, 2024
% Based on:"Multistep Nonlinear Predictive Controller
% David D. Brengel and Warren D. Seider, 1989.
% eq (28)-(29)
% Find the linearization arround the origin
syms x1 x2 u
xd = bioreact([x1 x2],u, para);
A = jacobian(xd, [x1 x2]);
A = eval(subs(A, {x1,x2,u}, {0,0,0}));
B = jacobian(xd, u);
B = eval(subs(B, {x1,x2,u}, {0,0,0}));
% gamma = para.gamma;
% CS1 = para.C1S; CS2 = para.C2S;
% beta= para.beta;
% Das = para.DaS;
% A = [                         - Das*exp(CS2/gamma)*(CS2 - 1) - 1,                                                                                                                                    - CS1*Das*exp(CS2/gamma) - (CS1*Das*exp(CS2/gamma)*(CS2 - 1))/gamma;
% -(Das*exp(CS2/gamma)*(CS2 - 1)*(beta + 1))/(beta - CS2 + 1), - (CS1*Das*exp(CS2/gamma)*(beta + 1))/(beta - CS2 + 1) - (CS1*Das*exp(CS2/gamma)*(CS2 - 1)*(beta + 1))/(beta - CS2 + 1)^2 - (CS1*Das*exp(CS2/gamma)*(CS2 - 1)*(beta + 1))/(gamma*(beta - CS2 + 1)) - 1];
% B = [-CS1
% -CS2];

end