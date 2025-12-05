function coeff = SolveCLF_new(A_in_eq,b_in_eq)
N  = size(A_in_eq,2);
c=ones(N,1);
coeff = linprog(c,A_in_eq,b_in_eq,...
    [],[],zeros(1,N),[],...
    optimoptions('linprog','Algorithm','interior-point','Display','off',...
    'MaxIterations',5000))';
% coeff = linprog(c,A_in_eq,b_in_eq)';
if numel(coeff) == 0
    error('infeasible LP')
end
end