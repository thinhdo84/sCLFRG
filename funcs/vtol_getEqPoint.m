function [xe, ue] = vtol_getEqPoint(x0,u0,epsilon)
[Res,~,exitflag,~]  =fsolve(@(X) vtol(X(1:6),X(7:8),epsilon),[x0;u0], ...
    optimoptions('fsolve','Display','off','algorithm','Levenberg-Marquardt'));
GreenFlags = [1 2 3 4];
if ismember(exitflag,GreenFlags)
    xe = round(Res(1:6),10);
    ue = round(Res(7:8),10);
    disp('Equilibrium found, vector field at equilibrium point:')
    disp(vtol(xe, ue, epsilon)');
    disp('----------------------------------------------------------------')
else
    xe = -999999;
    ue = 10000;
    warning('Equilibrium not found')
end
end