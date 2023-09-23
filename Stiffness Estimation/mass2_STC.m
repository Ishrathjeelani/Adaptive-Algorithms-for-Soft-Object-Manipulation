clear; clc;
close all;

global m1 m2 ts Ki stiff kp kd ts lambda l
ts=0.01; t0=0; tf=50;
tseries = t0:0.01:tf;
m1=2; m2=2; l=20; Ki=[5 0; 0 5]; kp=50; kd=5;

syms t
y1i=0; y1f=1; y2i=0; y2f=1;
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

yd=[y1f; y2f]; dyd=[dy1f; dy2f]; ddyd=[ddy1f; ddy2f]; stiff=2;

Fdes=[stiff*(y1f+y2f);stiff*(y1f+y2f)]; Fd=[rand(length(tseries),1)';rand(length(tseries),1)'];%Fd=[zeros(1,length(tseries)); zeros(1,length(tseries))];%
Fobs=[]; yact=[]; params=[0.1 0.1 0.1 0.1]'; Fact=[]; estimates=[];
lambda =0.01; Pest =inv(0.1*eye(4));

F=[0;0]; y=[0; 0]; dy=[0; 0]; P_hat=[0.001; 0.001]; Fe=[0; 0];errs=[]; Fest=[0; 0];estFor=[]; mest=[0.1 0.1]; kest=0.1;
for i=1:length(tseries)
    i
   [Fa,F,y,dy,ddy,e]=system(F,Fd(:,i),Fe,mest,kest,yd(:,i),dyd(:,i),ddyd(:,i),y,dy);
   [Fe,P_hat]=dob(P_hat,mest,kest,F,dy,y);
   Fe=-Fe;
   [Fest,params,forceErr,Pest]=estimateK(ddy,y,F-Fe,params,Pest);
   mest=[params(1) params(3)];
   kest=params(2);
   Fobs(:,end+1)=Fe;
   yact(:,end+1)=y;
   errs(:,end+1)=e;
   estFor(:,end+1)=Fest;
   Fact(:,end+1)=Fa;
   estimates(:,end+1)=params;
end
fig4=figure();
subplot(1,2,1);
plot(Fobs(1,:));
hold on;
plot(Fd(1,:))
legend('Estimated Disturbance','Actual Disturbance');
title("Disturbance Observer");
xlabel("Time(s)");
ylabel("Force(N)");
subplot(1,2,2);
plot(Fobs(2,:));
hold on;
plot(Fd(2,:))
legend('Estimated Disturbance','Actual Disturbance');
title("Disturbance Observer");
xlabel("Time(s)");
ylabel("Force(N)");

fig1=figure();
subplot(1,2,1);
plot(Fdes(1,:));
hold on;
plot(Fact(1,:));
legend('Force Desired','Actual Force');
title("Force on mass1 vs Time");
xlabel("Time(s)");
ylabel("Force(N)");

subplot(1,2,2);
plot(Fdes(2,:));
hold on;
plot(Fact(2,:));
legend('Force Desired','Actual Force');
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
subplot(1,3,1);
plot(estimates(1,:));
hold on;
plot(m1+zeros(1,length(tseries)));
legend('Estimated mass1','Real mass1');
title("Mass1 vs Time");
xlabel("Time(s)");
ylabel("Mass(Kg)");
subplot(1,3,2);
plot(estimates(3,:));
hold on;
plot(m2+zeros(1,length(tseries)));
legend('Estimated mass2','Real mass2');
title("Mass2 vs Time");
xlabel("Time(s)");
ylabel("Mass(Kg)");
subplot(1,3,3);
plot(estimates(4,:));
hold on;
plot(stiff+zeros(1,length(tseries)));
legend('Estimated Stiffness','Real Stiffness');
title("Stiffness vs Time");
xlabel("Time(s)");
ylabel("Stiffness(N/m)");

function [y_est,stiffi,e,P]=estimateK(ddx,x,y,stiffi,P)
global m1 m2 ts Ki stiff kp kd lambda l
X = [ddx(1) x(1)+x(2) ddx(2) x(1)+x(2)]';
K = (P*X) / (lambda + X'*P*X);
% K = (P*X)*inv((lambda + X'*P*X));
stiffi1 = [stiffi(1); stiffi(2)];
y_est1 = [ddx(1) x(1)+x(2)] * stiffi1;
e1 = y(1)-y_est1;
stiffi2 = [stiffi(3); stiffi(4)];
y_est2 = [ddx(2) x(1)+x(2)] * stiffi2;
e2 = y(2)-y_est2;
e=[e1; e1; e2; e2]';
y_est=[y_est1; y_est2];
stiffi_u1 =stiffi + K*e;
stiffi_u2 =stiffi + K*e;
stiffi = [stiffi_u1(1); stiffi_u1(2); stiffi_u2(1); stiffi_u2(2)];
P = P - K * X' * P;
end

function [Fa,F,y,dy,ddy,e]=system(F,Fd,Fe,mest,kest,yd,dyd,ddyd,y,dy)
global m1 m2 ts Ki stiff kp kd l
matM=[m1 0; 0 m2];
matK=[stiff stiff; stiff stiff];
ddy = inv(matM)*(F-matK*y+Fd);%-
%Fa=stiff*[(y(1)+y(2)); -(y(1)+y(2))]+Fd;
Fa=matK*[y(1); (y(2))]+Fd;%matM*ddy+
e=[y(1)-yd(1); y(2)-yd(2)];
de=dy-dyd;
dy=dy+ddy*ts;
y=y+dy*ts;
%y(2)=l-y(2);
F=[mest(1) 0; 0 mest(2)]*(ddyd-kp*e-kd*de)+kest*y-Fe;%;
end

function [Fe,P_hat]=dob(P_hat,mest,kest,F,dy,y)
global m1 m2 ts Ki stiff kp kd ts l
F=F-kest*(y(1)+y(2));
matM=[mest(1) 0; 0 mest(2)];
P = matM*dy;
err = P-P_hat;
P_hat= (F+Ki*err)*ts + P_hat;
Fe = Ki*(P-P_hat);
%Fe = Fe+kest*(y(1)+y(2));
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