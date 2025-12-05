function alph_result=Ellipse_inner_fixed_shape(P,Qp)
% Find the maximum volume ellipsoid E inside a convex polytopic set P
% given a fixed matrix Qp as in x'*Qp*x<=alpha
% 
% Details can be found here
% https://web.stanford.edu/~boyd/cvxbook/bv_cvxslides.pdf  % Slide 184
% https://math.stackexchange.com/questions/509486/biggest-ellipse-included-in-a-convex-polygon
% https://see.stanford.edu/materials/lsoeldsee263/15-symm.pdf %SVD
% Input: the bounding polytope P
% Output: alpha
%
% Developer:
% ----------------------------------------------------------------------------------
% Huu Thinh DO

% n=size(P.A,2);
sdpvar sqr_alp bet;
[U,S,V] = svd(Qp);
J=U';
Py=Polyhedron(P.A*J',P.b);
Py=P;
constr=sqr_alp>=0;
% MM = (V * (S)^(1/2) *U')^-1;
MM = ((U*sqrt(S)*V')')^-1;
B = sqr_alp * MM;
% obj=log(det(B));
obj = logdet(B);
% obj = 0;
for i=1:size(Py.A,1)
    constr=[constr, (B*Py.A(i,:)')'*(B*Py.A(i,:)')  <=  (Py.b(i,:))^2];
end
options=sdpsettings('solver','sdpt3','sdpt3.maxit',4000,'verbose',0);
diagnostics = optimize(constr,-obj,options);
%% Result
if diagnostics.problem==0
    disp(diagnostics)
    disp('Largest ellipsoid found')
    alph_result=value(sqr_alp)^2;
else 
    disp(diagnostics)
    error 'cannot find the solution'
end
end
