function plot_vtol(x,y,theta,d,colorPara)
alph = colorPara.patchalpha;
col = colorPara.patchcolor;
lc = colorPara.linecolor;
Rot = [cos(theta), -sin(theta); sin(theta), cos(theta)];
Rot_small = [cos(theta), sin(theta); -sin(theta), cos(theta)];

r = 0.7*d/3;

Pleft = [0, r;
        0, 0.5*r;
        -0.85*d, -0.32*r;
        -d,0];
Pright = [0, r;
         0, 0.5*r;
        0.85*d, -0.32*r;
        d,0];
Cir = random_points_on_ellipsoid([1, 0; 0, 1.5],(0.5*r)^2,[0;0.8*r],12);
for i =1:size(Cir,2)
Cir(:,i) = Rot*Cir(:,i) + [x;y];
end
for i =1:size(Pleft,1)
Pleft(i,:) = Pleft(i,:)*Rot'  + [x,y];
Pright(i,:) = Pright(i,:)*Rot' + [x,y];
end
patch(Pleft(:,1), Pleft(:,2),col,'facealpha',alph,'edgecolor',lc);
hold on
patch(Pright(:,1), Pright(:,2),col,'facealpha',alph,'edgecolor',lc);
patch(Cir(1,:), Cir(2,:),col,'facealpha',alph,'edgecolor',lc);
Ellipsoid_patch2D(Rot*[1.5, 0; 0, 1],(r)^2,[x;y],'edgecolor',lc,...
    'facecolor',col,'facealpha',alph);



end