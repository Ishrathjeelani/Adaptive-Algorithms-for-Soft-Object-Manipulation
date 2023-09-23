clear all; close all; clc;

k0=2; k=1; theta=1; r=0.1;

ts=0.1; ti=0; te=50; 
t = ti:ts:te;
uc = 5*sin(t);
Gm=tf(k0,[1 1]);
ym=lsim(Gm,uc,t);
res=[];

for i=1:length(t)
G=tf(k,[1 1]);
u=theta.*uc';
y=lsim(G,u,t);
e=y-ym;
A=tf(-r,[1 0]);
theta=lsim(A,ym.*e,t);
%theta=theta(i);
u=theta.*uc';
y=lsim(G,u,t);
res(end+1)=y(i);
end


plot(ym);
hold on;
plot(res);
legend('Model Response','Actual Response');
title("MRAC Input/Output signals");
xlabel("Time(s)");
ylabel("Amplitude");