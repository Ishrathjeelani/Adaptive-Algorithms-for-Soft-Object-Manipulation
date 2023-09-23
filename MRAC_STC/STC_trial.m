% clear; clc;
% close all;

global m1 ts Ki stiff kp kd ts lambda
ts=0.01; t0=0; tf=20;
tseries = t0:0.01:tf;
m1=5; l=20; kp=5; kd=0.1; stiff=2; lambda =0.01; Pest =inv(lambda*eye(2)); Ki=5;

syms t
y1i=0; y1ff=2; 
[y1d,dy1d,ddy1d]= trajectory_planner(t0,tf,y1i,0,y1ff,0);
for i=1:1:length(tseries)
    y1f(i) = double(subs(y1d,t,tseries(i)));
    dy1f(i) = double(subs(dy1d,t,tseries(i)));
    ddy1f(i) = double(subs(ddy1d,t,tseries(i)));
end
Fdes=m1*ddy1f+stiff*y1f; Fd=0.1+zeros(length(tseries));

kest=0; mest=1; e=0; de=0; P=Pest; params=[m1; kest]; y=y1i; dy=0; F_des=0; Fe=0; P_hat=0;
actPath=[]; mTseries=[]; kTseries=[];forces=[]; act_force=[]; obsForce=[];
for i=1:length(tseries)
[F]=controller(m1,kest,Fe,ddy1f(i),y1f(i),e,de);
[e,de,y,dy,ddy,Fext]=system(F,Fd(i),y1f(i),dy1f(i),ddy1f(i),y,dy,params);
[Fe,P_hat]=dob(P_hat,F,dy);
[y_est, P, err, params]=estimateParams(y,ddy,(m1*ddy1f(i)+Fe),y1f(i),ddy1f(i),P,params);%ideally we don't know Fdes
mest=params(1);
kest=params(2);
params=[m1; kest];
actPath(end+1)=y;
mTseries(end+1)=mest;
kTseries(end+1)=kest;
forces(end+1)=-y_est;
act_force(end+1)=Fext;
obsForce(end+1)=Fe;
end

% fig4=figure();
% plot(obsForce);
% hold on; 
% plot(-stiff*(y1f));

fig1= figure();
plot(tseries,actPath,'b');
hold on;
plot(tseries,y1f,'r');
legend('Desired Trajectory','Actual Trajectory');
title("Trajectory Tracking");
xlabel("Time(s)");
ylabel("Distance(mm)");
hold off;


fig2= figure();
subplot(2,1,1);
hold on;
plot(tseries,mTseries,'b');
plot(tseries,m1+zeros(length(tseries)),'r');
legend('Estimated Mass','Actual Mass');
title("Mass Estimation");
xlabel("Time(s)");
ylabel("Mass(Kg)");
hold off;
subplot(2,1,2);
plot(tseries,-kTseries,'b');
hold on;
plot(tseries,stiff+zeros(length(tseries)),'r');
legend('Estimated Stiffness','Actual Stiffness');
title("Stiffness Estimation");
xlabel("Time(s)");
ylabel("Stiffness(N/m)");
hold off;

fig3= figure();
plot(tseries,forces,'b');
hold on;
plot(tseries,Fdes,'r');
plot(tseries,act_force,'g');
legend('Estimated Forces','Desired Forces','Actual Forces');
title("Force Tracking");
xlabel("Time(s)");
ylabel("Force(N)");
hold off;

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

function [y_est, P, err,params]=estimateParams(x,ddx,y,y1f,ddy1f,P,params)
global m1 m2 ts Ki stiff kp kd lambda
X = [ddx x]';
K = (P*X) / (lambda + X'*P*X);
% K = (P*X)*inv((lambda + X'*P*X));
y_est = X' * params;
%y = params(2)*ddy1f-params(4)*y1f;
err = y-y_est;
params =params + K*err;
P = P - K * X' * P;
end

function [F]=controller(mest,kest,Fe,ddyd,yd,e,de)
global m1 ts Ki stiff kp kd
F=mest*(ddyd-kp*e-kd*de)-Fe;%;+kest*yd
end

function [Fe,P_hat]=dob(P_hat,F,dy)
global m1 ts Ki stiff kp kd
matM=m1;
Pot = matM*dy;
err = Pot-P_hat;
P_hat= (F+Ki*err)*ts + P_hat;
Fe = Ki*(Pot-P_hat);
end

function [e,de,y,dy,ddy,Fext]=system(F,Fd,yd,dyd,ddyd,y,dy,param)
global m1 ts Ki stiff kp kd
% estM=param(1);
% estK=param(2);
ddy = (F+Fd-stiff*y)/m1; %(F+Fd+stiff*y)/m1;%-
e=y-yd;
de=dy-dyd;
dy=dy+ddy*ts;
y=y+dy*ts;
Fext=m1*ddy+stiff*y;
end