clc
clear
close all
%%
tau = 0.01;
tsim = 10;
Nsim = round(tsim/tau);
Xsim = zeros(2,Nsim+1);
usim = zeros(1,Nsim) + 97;
Xsim(:,1) = [0.0795; 443.4566];
for i =1:Nsim
    Xsim(:,i+1) = Xsim(:,i) + tau*cstr(Xsim(:,i),usim(:,i));
end

%%
figure
subplot(2,1,1)
plot(Xsim(1,:))
subplot(2,1,2)
plot(Xsim(2,:))

