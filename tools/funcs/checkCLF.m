function res = checkCLF(V,A,B,beta,U,X)
% check the CLF conditions
% V is the parameteried Lyapunov object
% U is the constraint set
% the wrong version
% check = V.grad(X) * (A*X  +  B*U.support(-B'*V.grad(X)') ) + beta * V.eval(X);
% the correct one?
check = V.grad(X) * A*X  - U.support(-B'*V.grad(X)')  + beta * V.eval(X);
% U.support(-B'*V.grad(X)')
% U.support(B'*V.grad(X)')
% if U.support(-B'*V.grad(X)') < -1e-8
%     disp(U.support(-B'*V.grad(X)'))
%     error('this value is not supposed to be nagative')
% end
res = check <=0;
% check
end