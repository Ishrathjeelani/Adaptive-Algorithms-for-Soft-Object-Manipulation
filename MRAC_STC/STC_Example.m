clear all; close all; clc;
% True System Params
a=5; b=1;

% Reference System Params
w=4; zeta=0.8;


G_s=tf([b],[1 a 0]);
Gm_s=tf([w^2],[1 2*zeta*w w^2]);

ts=0.1; ti=0; te=20; b0=1;
t = ti:ts:te;
uc = square(t);%b0*sin(2*t);
y0=0;
ym=lsim(Gm_s,uc,t);
y=lsim(G_s,uc,t);

t0=0; s0=0; a0=2; s1=0; r1=0; u=uc; dy0=0; ddy0=1;
res=[]; estimates=[]; subset=1; P=100*eye(2); th_e=[5; 10];

for i=1:length(t)-subset
pltf=tf([b],[1 a 0]);
yc=lsim(pltf,u,t);
fftf = tf([t0 t0*a0],[1 r1]);
u1=lsim(fftf,uc,t);
fbtf = tf([s0 s1],[1 r1]);
u2=lsim(fbtf,yc,t);
u=u1-u2;

dy=gradient(yc,t);
ddy=gradient(dy,t);
if(rem(i,subset)==0)
    phi=[]; Y=[]; K=[];
    for j=i:i+subset-1
        phi(end+1,:)=[ddy(j) dy(j)];
        Y(end+1)=u(j);
    end
K=P*phi'*inv(eye(length(phi*P*phi'))+phi*P*phi');%P*phi';
err=Y'-phi*th_e;
P=(eye(2)-K*phi)*P;
% phi=[1 u(i+1) u(i+1)^2; 1 u(i+2) u(i+2)^2; 1 u(i+3) u(i+3)^2];
% Y=[ym(i+1) ym(i+2) ym(i+3)]';
th_e = th_e+K*(Y'-phi*th_e);
%th_e=linsolve(phi,Y');
%th_e=inv(phi)*Y';
a_h=th_e(2)/th_e(1);%param(2)/param(1);
b_h=1/th_e(1);%1/param(2);
estimates(:,end+1)=[a_h; b_h];
r1=2*zeta*w+a0-a;
s0=(a0*2*zeta*w+w^2-a*r1)/b;
s1=(w^2*a0)/b;
t0=w^2/b;
end
%[P_mat,param,r1,s0,s1,t0] = estimator(u(i),param,P_mat,lambda,zeta,w,a0,dy(i),ddy(i));
res(end+1)=yc(i);
end

fig1=figure();
plot(uc);
hold on;
plot(res);
plot(ym);
%plot(y);
legend('Input Signal','Actual Response', 'Model Response');
hold off;
title("STC Input/Output signals");
xlabel("Time(s)");
ylabel("Amplitude");

fig2=figure();
subplot(2,1,1);
plot(a+zeros(length(t)));
hold on;
plot(estimates(1,:));
legend('Actual Value','Estimated Value');
hold off;
title("Estimation of Parameter a");
xlabel("Time(s)");
ylabel("Value");

subplot(2,1,2);
plot(b+zeros(length(t)));
hold on;
plot(estimates(2,:));
legend('Actual Value','Estimated Value');
hold off;
title("Estimation of Parameter b");
xlabel("Time(s)");
ylabel("Value");

