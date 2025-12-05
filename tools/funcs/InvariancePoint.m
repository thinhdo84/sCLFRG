function res = InvariancePoint(V,fx,gx,U,X)
% check the invariant condition at the point X
check = V.grad(X)*fx-U.support(-gx'*V.grad(X)') ;
res = check <=0;
end