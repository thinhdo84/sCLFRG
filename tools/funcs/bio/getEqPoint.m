function CS = getEqPoint(x0,para)
[cs_tmp,~,exitflag,~]  = fsolve(@(Ceq) massBalance(Ceq, para.beta, para.gamma, para.DaS),x0,...
    optimoptions('fsolve','Display','off'));
GreenFlags = [1 2 3 4];
if ismember(exitflag,GreenFlags)
    CS = cs_tmp;
    disp('Equilibrium point found')
    disp('The state derivetive:')
    disp(massBalance(CS, para.beta, para.gamma, para.DaS))
else
    CS = 9999 + 8888i;
    error('Equilibrium point not found')
end
end

function Cd = massBalance(x, beta, gamma, DaS)

Cd(1) = -x(1) + DaS* (1-x(2))*exp(x(2)/gamma)*x(1);
Cd(2) = -x(2) + DaS*(1+beta)*(1-x(2))*exp(x(2)/gamma)*x(1)/(1 + beta - x(2));
end


