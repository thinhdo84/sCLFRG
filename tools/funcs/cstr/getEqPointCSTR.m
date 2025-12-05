function y = getEqPointCSTR(u, x_guess)
y=fsolve(@(xeq) cstr(xeq,u),x_guess, optimoptions('fsolve','Display','off'));
fprintf('Equilibrium point found for input u = %.3f\n',u)
end