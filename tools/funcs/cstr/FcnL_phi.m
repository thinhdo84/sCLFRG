% function L_Phi = FcnL_phi(AK,K,P,alpha,xe,ue)
%     opt = optimset('MaxFunEvals',10000,'MaxIter',10000,'Display','off');
%     objF = @(xt) -sqrt(FcnPhi(xt,xe,ue,AK,K)' * FcnPhi(xt,xe,ue,AK,K))/sqrt(xt'*xt);
%     [~,L_Phi_tilde] = fmincon( objF,...
%         xe,[],[],[],[],[],[],@(xt)nonlinConsAlpha(xt,P,alpha),opt);
%     L_Phi = -L_Phi_tilde;
% end
function L_Phi = FcnL_phi(AK,K,P,alpha,xe,ue)
pts = random_points_on_ellipsoid(P,alpha,[0;0],20);
ell_approx = Polyhedron('V',pts');
A = ell_approx.A; b=ell_approx.b;
opt = optimset('MaxFunEvals',10000,'MaxIter',10000,'Display','off');

[~,L_Phi_tilde] = fmincon(  @(xt) -sqrt(FcnPhi(xt,xe,ue,AK,K)' * FcnPhi(xt,xe,ue,AK,K))/sqrt(xt'*xt),...
        xe,A,b,[],[],[],[],[],opt);
L_Phi = -L_Phi_tilde;
end