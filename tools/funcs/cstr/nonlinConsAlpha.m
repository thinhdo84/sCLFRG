function [c, ceq] = nonlinConsAlpha(xt, P, alpha)
    c = xt'*P*xt - alpha;
    ceq = [];
end
