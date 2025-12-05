function u = solveMPC(x,xinit,ut,solver)
try
    solver.set_value(xinit,x);
    sol=solver.solve();
    u=sol.value(ut(:,1));
catch
    m = size(ut,1);
    u = nan(m,1);
end
end