function p = FcnPhi(x,AK,K,para)
% difference between the jacobian linearization and the nonlinear function
p = bioreact(x,K*x, para) - AK*x;
end