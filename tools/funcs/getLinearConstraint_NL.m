function [ai, bi]=getLinearConstraint_NL(V,dynamics,xi,ui,g)
N = numel(V); 
n = numel(xi);
Vtmp = zeros(n,N);
for i =1:N
    Vtmp(:,i) = V{i}.grad(xi)'; 
end
ai = dynamics(xi,ui)'*Vtmp;
bi = -g;
% ai = zeros(1,N);
% for i =1:N
%     ai(i) = V{i}.grad(xi)*dynamics(xi,ui); 
% end
end
