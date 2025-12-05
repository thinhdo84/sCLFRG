function u_sat = sat_polytope(u,U_constraint)
% Saturation function for an arbitrary polytopic set
if all(all(U_constraint.A*u <= U_constraint.b))
    u_sat=u;
else
    lambda_max = 1e5; %some large number
    for i =1:size(U_constraint.A,1)
        if U_constraint.A(i,:)*u>0
            lambda_max = min(lambda_max,min(U_constraint.b(i)/(U_constraint.A(i,:)*u),1));
        end
    end
    u_sat = u * lambda_max;
end
end