function [coeff, msg] = SolveCLF(Vbasis, weight, gamma, gridX, gridU, g_cost, sys)
N = numel(Vbasis);
alph = sdpvar(1,N);
M = size(gridX,2); % the number of pairs
% Q=g_cost.Q; R = g_cost.R;
sumalpha=0;
objective = 0;
constraints = [alph>=0];
for k = 1:N
    objective = objective + weight(k)* alph(k);
end
% objective =alph*alph';
for i = 1:M
    lya_cond_tmp = 0;
    X_i = gridX(:,i);
    u_i= gridU(:,i);
    for k = 1:N
        lya_cond_tmp = lya_cond_tmp +...
            alph(k) *(Vbasis{k}.grad(X_i)*(sys.A*X_i+ sys.B*u_i)) +...
            gamma*g_cost(i);
    end
%     lya_cond_tmp
%     i
    constraints = [constraints, lya_cond_tmp <= 0]; 
end

options=sdpsettings('solver','linprog');
disp('constraints formulated')
msg = optimize(constraints,objective,options);
if msg.problem >0
    disp(msg)
    error('sub CLF not found')
else
    disp('sub-CLF candidate found')
end
coeff=value(alph);



end