clear all; close all; clc;

ts=0.1; ti=0; te=20; 
t = ti:ts:te;
uc = 5*sin(t);
z = uc + 0.8*rand(length(t),1)';

x=[0]; P=15; Q=1; R=10; H=1;
for i=2:length(t)
xp = x(i-1);
P = P + Q;
K = P*H'*inv(H*P*H'+ R);
x(i) = xp + K*(z(i)-H*xp);
P = P - K*H*P;
end

plot(z);
hold on;
plot(x);