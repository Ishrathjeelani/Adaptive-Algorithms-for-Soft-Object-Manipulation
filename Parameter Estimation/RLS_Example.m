close all; clear all; clc;

%% 
% Least square estimation with Gaussian Noise/White Noise.
t = (0:0.1:60)';
u = sin(t);%sawtooth(t);

b0=1.2; b1=1.3; b2=1.5;
y1= b0+b1*u+b2*u.^2; std=0.1; meanValue=0;
ym = b0+b1*u+b2*u.^2;%+ std*randn(size(y1)) + meanValue;%+awgn(y1,10,'measured');
estimates=[]; subset=1; P=100*eye(3); th_e=[0.1; 0.1; 0.1];
for i=1:subset:length(t)-subset
    i
    phi=[]; Y=[]; K=[];
    for j=i:i+subset-1
        phi(end+1,:)=[1 u(j) u(j)^2];
        Y(end+1)=ym(j);
    end
K=P*phi'*inv(eye(length(phi*P*phi'))+phi*P*phi');%P*phi';
err=Y'-phi*th_e;
P=(eye(3)-K*phi)*P;
% phi=[1 u(i+1) u(i+1)^2; 1 u(i+2) u(i+2)^2; 1 u(i+3) u(i+3)^2];
% Y=[ym(i+1) ym(i+2) ym(i+3)]';
th_e = th_e+K*(Y'-phi*th_e);
%th_e=linsolve(phi,Y');
%th_e=inv(phi)*Y';
estimates(:,end+1)=th_e;
end
l=length(estimates);
fig1=figure();
plot(ym);
fig2=figure();
subplot(3,1,1);
plot(b0+zeros(l));
hold on; 
plot(estimates(1,:));
subplot(3,1,2);
plot(b1+zeros(l));
hold on; 
plot(estimates(2,:));
subplot(3,1,3);
plot(b2+zeros(l));
hold on; 
plot(estimates(3,:));
