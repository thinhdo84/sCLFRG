function [cieq, ceq ] = NLconstraint(x,P,a)
cieq=x'*P*x - a; ceq = [];
end