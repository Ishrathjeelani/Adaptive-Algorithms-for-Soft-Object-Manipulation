clear; clc;
close all;

global m1 m2 ts Ki stiff kp kd ts lambda
ts=0.01; t0=0; tf=5;
tseries = t0:0.01:tf;
m1=2; m2=4;l=20; Ki=[12 0; 0 12]; kp=100; kd=0.1;

syms t
y1i=0; y1f=2; y2i=0; y2f=-3;
[y1d,dy1d,ddy1d]= trajectory_planner(t0,tf,y1i,0,y1f,0);
[y2d,dy2d,ddy2d]= trajectory_planner(t0,tf,y2i,0,y2f,0);
for i=1:1:length(tseries)
    y1f(i) = double(subs(y1d,t,tseries(i)));
    y2f(i) = double(subs(y2d,t,tseries(i)));
    dy1f(i) = double(subs(dy1d,t,tseries(i)));
    dy2f(i) = double(subs(dy2d,t,tseries(i)));
    ddy1f(i) = double(subs(ddy1d,t,tseries(i)));
    ddy2f(i) = double(subs(ddy2d,t,tseries(i)));
end

% fig3=figure();
% plot(-stiff*(y1f-y2f));

yd=[y1f; y2f]; dyd=[dy1f; dy2f]; ddyd=[ddy1f; ddy2f]; stiff=6*(y1f-y2f);

Fd=[stiff.*(y1f-y2f);-stiff.*(y1f-y2f)];%Fd=[zeros(1,length(tseries)); zeros(1,length(tseries))];%
Fes=[]; yact=[]; stiff_est=[0.1 0.1 0.1 0.1]';
lambda =0.1; Pest =inv(lambda*eye(4));

F=[0;0]; y=[0; 0]; dy=[0; 0]; P_hat=[0; 0]; Fe=[0; 0];errs=[]; Fest=[0; 0];Kest=[];
for i=1:length(tseries)
   [F,y,dy,e]=system(F,Fd(:,i),Fe,stiff_est,yd(:,i),dyd(:,i),ddyd(:,i),y,dy);
   [Fe,P_hat]=dob(P_hat,F,dy);
   [Fest,stiff_est,forceErr,Pest]=estimateK(y(1)-y(2),Fd(1,i),stiff_est,Pest);
   Fes(:,end+1)=Fe;
   yact(:,end+1)=y;
   errs(:,end+1)=e;
   Kest(:,end+1)=Fest;
end

fig1=figure();
subplot(1,2,1);
plot(Fes(1,:));
hold on;
plot(Fd(1,:));
legend('Force Estimated','Actual Force');
title("Force on mass1 vs Time");
xlabel("Time(s)");
ylabel("Force(N)");

subplot(1,2,2);
plot(Fes(2,:));
hold on;
plot(Fd(2,:));
legend('Force Estimated','Actual Force');
title("Force on mass2 vs Time");
xlabel("Time(s)");
ylabel("Force(N)");

fig2=figure();
subplot(1,2,1);
plot(yact(1,:));
hold on;
plot(yd(1,:));
legend('Actual Trajectory','Desired Trajectory');
title("Trajectory of mass1 vs Time");
xlabel("Time(s)");
ylabel("Position(m)");

subplot(1,2,2);
plot(yact(2,:));
hold on;
plot(yd(2,:));
legend('Actual Trajectory','Desired Trajectory');
title("Trajectory of mass2 vs Time");
xlabel("Time(s)");
ylabel("Position(m)");
%y_predicted = [y_predicted; y_est];

fig3=figure();
plot(Kest(1,:)./(y1f-y2f));
hold on;
plot(stiff);
legend('Estimated Stiffness','Real Stiffness');
title("Stiffness vs Time");
xlabel("Time(s)");
ylabel("Stiffness(N/m)");

function [y_est,stiffi,e,P]=estimateK(x,y,stiffi,P)
global m1 m2 ts Ki stiff kp kd lambda
X = [1 x x^2 x^3]';
K = (P*X) / (lambda + X'*P*X);
% K = (P*X)*inv((lambda + X'*P*X));
y_est = X' * stiffi;
e = y-y_est;
stiffi =stiffi + K*e;
P = P - K * X' * P;
end

function [F,y,dy,e]=system(F,Fd,Fe,stiff_est,yd,dyd,ddyd,y,dy)
global m1 m2 ts Ki stiff kp kd
matM=[m1 0; 0 m2];
matK=[stiff_est(1) stiff_est(1); stiff_est(2) stiff_est(2)];
ddy = inv(matM)*(F+Fd);%-
e=y-yd;
de=dy-dyd;
dy=dy+ddy*ts;
y=y+dy*ts;
F=matM*(ddyd-kp*e-kd*de)-Fe;%;
end

function [Fe,P_hat]=dob(P_hat,F,dy)
global m1 m2 ts Ki stiff kp kd ts
matM=[m1 0; 0 m2];
P = matM*dy;
err = P-P_hat;
P_hat= (F+Ki*err)*ts + P_hat;
Fe = Ki*(P-P_hat);
end

function [yd,dyd,ddyd]= trajectory_planner(t0,tf,yi,dyi,yf,dyf)
syms t
Q = [yi;dyi;yf;dyf];
B = [1 t0 t0*t0 t0*t0*t0;
    0 1 2*t0 3*t0*t0;
    1 tf tf*tf tf*tf*tf;
    0 1 2*tf 3*tf*tf];
A1= inv(B)*Q;
a0=A1(1,1);
a1=A1(2,1);
a2=A1(3,1);
a3=A1(4,1);
yd=a0+(a1*t)+(a2*t*t)+(a3*t*t*t); % Thita(t)
dyd=(a1)+(2*a2*t)+(3*a3*t*t);
ddyd=2*a2+ (6*a3*t);
end