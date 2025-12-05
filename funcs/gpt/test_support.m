clc
clear
close all
%%


H = [0.3192    1.0933
    0.3129    1.1093
   -0.8649   -0.8637
   -0.0301    0.0774
   -0.1649   -1.2141
    0.6277   -1.1135];
p = 2;               % use 2p-norm with p = 2 (i.e., 4-norm)
% V = Lya
ComputeContour2D(rBound,V,levelValue)

sigma = support_functional_smoothPolytope(H, p, y);
disp(sigma);
