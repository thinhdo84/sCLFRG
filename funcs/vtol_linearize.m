function [A, B] = vtol_linearize(xe,ue,epsilon)
% Created by Huu-Thinh DO, 2024
syms x1 x2 x3 x4 x5 x6 u1 u2
X=[x1; x2; x3 ;x4 ;x5 ;x6];
u = [u1;u2];
eqpts = reshape([num2cell(xe); num2cell(ue)],[1,8]);
xd =  vtol(X,u,epsilon);
A = jacobian(xd,X);
A =eval( subs(A, {x1 x2 x3 x4 x5 x6 u1 u2}, eqpts));
B = jacobian(xd,u);
B = eval(subs(B, {x1 x2 x3 x4 x5 x6 u1 u2}, eqpts));
end

