function u = solveMPC_tracking(x,ref,xinit,xref,ut,solver)
try
    solver.set_value(xinit,x);
    solver.set_value(xref,ref);
    sol=solver.solve();
    u=sol.value(ut(:,1));
catch
    m = size(ut,1);
    u = nan(m,1);
end
end