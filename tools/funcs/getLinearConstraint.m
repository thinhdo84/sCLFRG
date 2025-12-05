function [ai, bi]=getLinearConstraint(V,A,B,xi,ui,g)
N = numel(V); n = size(B,1);
Vtmp = zeros(n,N);
for i =1:N
    Vtmp(:,i) = V{i}.grad(xi)'; 
end
ai = (A*xi + B*ui)'*Vtmp;
bi = -g;
end
