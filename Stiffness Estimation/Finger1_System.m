clear; clc;
close all;

global m1 m2 ts Ki kp kd l1 l2 l0 lc1 lc2 I1 I2 g stiff ki lambda
ts=0.01; t0=0; tf=5;
tseries = t0:ts:tf;
m1=2; m2=2; l0=20; Ki=[2.5 0; 0 2.5]; kp=50; kd=5; ki=0; stiff=5; lambda =0.01;

%%Finger Dimension
l1=1.5; l2=1.5; lc1=l1/2; lc2=l2/2; I1=l1^2/3; I2=l2^2/3; g=9.81;

syms t
x1i=0; x1f=1; y1f=2;
[xd,dxd,ddxd]= trajectory_planner(t0,tf,x1i,0,x1f,0,y1f);
q1=[]; dq1=[]; q2=[]; dq2=[]; prevq1=0; prevq2=0; ddqd=[0;0];
%[y2d,dy2d,ddy2d]= trajectory_planner(t0,tf,y2i,0,y2f,0,z);
for i=1:1:length(tseries)
    i
    x1f(i) = double(subs(xd,t,tseries(i)));
    dx1f(i) = double(subs(dxd,t,tseries(i)));
    [x,dx]=inverse_kine(x1f(i),dx1f(i),y1f);
    ddqd(:,end+1)=[dx(1)-prevq1; dx(2)-prevq2]./ts;
    q1(end+1)=x(1); 
    dq1(end+1)=dx(1);
    prevq1= dx(1);
    q2(end+1)=x(2);
    dq2(end+1)=dx(2);
    prevq2= dx(2);
end
q=[q1;q2]; dq=[dq1; dq2];


T=[0;0]; Fd=[zeros(1,length(tseries));zeros(1,length(tseries))]; Te=[0;0]; y=q(:,1); dy=dq(:,1);%Fext=[stiff*y1f;zeros(1,length(tseries))];
P_hat=[0;0]; estF=[]; acty=[]; err=[0;0]; torq=[];
stiff_est=[0]'; xf=[]; yf=[]; xi=0; yi=0;
 Pest =inv(lambda*eye(1)); Kf=[];
for i=1:length(tseries)
 [T,y,dy,Text,err]=system(T,Fd(i),Te,stiff*[xi;yi],q(:,i),dq(:,i),ddqd(:,i),y,dy,err);
 [Fe,Te,P_hat]=dob(P_hat,T,y,dy);
 [xi,yi]=forward_kin(y);
 [y_est,stiff_est,e,Pest]=estimateK(xi,Fe,stiff_est,Pest);
 estF(:,end+1)=Fe;
 torq(:,end+1)=stiff*[xi;yi-y1f];%Fext(:,i);
 acty(:,end+1)=y;
 Kf(:,end+1)=stiff_est;
 xf(end+1)=xi;
 yf(end+1)=yi;
end


fig1=figure();
subplot(2,1,1);plot(q1);hold on;plot(acty(1,:));
legend('Desired Q1','Actual Q1');
title("Theta1 Trajectory Tracking");
xlabel("Time(s)");
ylabel("Angle(rad)");
subplot(2,1,2);plot(q2);hold on;plot(acty(2,:));
legend('Desired Q2','Actual Q2');
title("Theta2 Trajectory Tracking");
xlabel("Time(s)");
ylabel("Angle(rad)");

fig4=figure();
subplot(2,1,1); plot(xf); hold on; plot(x1f);
legend('Actual Displacement','Desired Displacement');
title("Trajectory Tracking along x-axis");
xlabel("Time(s)");
ylabel("Displacement(mm)");
subplot(2,1,2); plot(yf); hold on; plot(y1f+zeros(length(tseries))); %
legend('Actual Displacement','Desired Displacement');
title("Trajectory Tracking along y-axis");
xlabel("Time(s)");
ylabel("Displacement(mm)");

fig2=figure();
subplot(2,1,1);plot(torq(1,:));hold on;plot(estF(1,:));
legend('Actual Force','Estimated Force');
title("Force Estimation along x-axis");
xlabel("Time(s)");
ylabel("Force(N)");
subplot(2,1,2);plot(torq(2,:));hold on;plot(estF(2,:));
legend('Actual Force','Estimated Force');
title("Force Estimation along y-axis");
xlabel("Time(s)");
ylabel("Force(N)");

fig3=figure();
plot(x1f,Kf(1,:)); hold on; plot(x1f,stiff+zeros(length(tseries)));
legend('Estimated Stiffness','Actual Stiffness');
title("Stiffness Estimation");
xlabel("Displacement(mm)");
ylabel("Stiffness(N/mm)");

function [x,z]=forward_kin(y)
global m1 m2 ts kp kd l1 l2 l0 lc1 lc2 I1 I2 g ki
q1i=y(1); q2i=y(2); 
x = l1*cos(q1i) + l2*cos(q2i);
z = +l1*sin(q1i) + l2*sin(q2i);
end

function [y_est,stiffi,e,P]=estimateK(x,y,stiffi,P)
global m1 m2 ts kp kd l1 l2 l0 lc1 lc2 lambda
X = [x]';
K = (P*X) / (lambda + X'*P*X);
% K = (P*X)*inv((lambda + X'*P*X));
y_est = X' * stiffi;
e = y-y_est;
stiffi =stiffi + K*e;
P = P - K * X' * P;
end

function [T,y,dy,Text,err]=system(T,Fd,Te,Fext,yd,dyd,ddqd,y,dy,err)
global m1 m2 ts kp kd l1 l2 l0 lc1 lc2 I1 I2 g ki
q1=y(1); q2=y(2); dq1=dy(1); dq2=dy(2);
% m11=m1*lc1^2+m2*(l1^2+lc1^2+2*l1*lc2*cos(q2))+I1+I2;
% m12=m2*(lc2^2+l1*lc2*cos(q2))+I2;
% m21=m2*(lc2^2+l1*lc2*cos(q2))+I2;
% m22=m2*lc2^2+I2;
% m11=m1*lc1^2+m2*l1^2+I1;
% m12=m2*l1*lc2*cos(q2);
% m21=m12;
% m22=m2*lc2^2+I2;
% M=[m11 m12;
%     m21 m22];
% C=-m2*l1*lc2*sin(q2)*g;%-m2*l1*lc2*sin(q2)*[2*dq1*dq2+dq2^2; -dq1^2];

%% Dynamics
M = [(m1*l1^2/4)+m2*l1^2+I1, (m2*l1*l2*cos(q1-q2))/2;
          (m2*l1*l2*cos(q1-q2))/2, (m2*l2^2)/4 + I2] ;  
C = [0,(m2*l1*l2*sin(q2-q1))/2;
     -(m2*l1*l2*sin(q2-q1))/2, 0] ;
G=0; %[(m1*lc1+m2*l1)*g*cos(q1); m2*g*lc2*cos(q1+q2)];%[(m1*lc1+m2*l1)*g*cos(q1)+m2*g*lc2*cos(q1+q2); m2*g*lc2*cos(q1+q2)];
J=[-l1*sin(q1) -l2*sin(q2); l1*cos(q1) l2*cos(q2)];
Text=transpose(J)*Fext;
ddy=inv(M)*(T+Text-C*dy-G);
e=y-yd;
err=err+e;
de=dy-dyd;
dy=dy+ddy*ts;
y=y+dy*ts;
T=M*(ddqd-kp*e-kd*de-ki*err)+C*dy-Te;%;ddyd
end

function [Fe,Te,P_hat]=dob(P_hat,T,y,dy)
global m1 m2 ts kp kd l1 l2 l0 lc1 lc2 I1 I2 g Ki
q1=y(1); q2=y(2); dq1=dy(1); dq2=dy(2);
% m11=m1*lc1^2+m2*l1^2+I1;
% m12=m2*l1*lc2*cos(q2);
% m21=m12;
% m22=m2*lc2^2+I2;
% matM=[m11 m12;
%     m21 m22];
% C=-m2*l1*lc2*sin(q2)*g;
matM = [(m1*l1^2/4)+m2*l1^2+I1, (m2*l1*l2*cos(q1-q2))/2;
          (m2*l1*l2*cos(q1-q2))/2, (m2*l2^2)/4 + I2] ;   
C = [0,(m2*l1*l2*sin(q2-q1))/2;
     -(m2*l1*l2*sin(q2-q1))/2, 0] ;
Pot = matM*dy;
err = Pot-P_hat;
P_hat= (T+C'*dy+Ki*err)*ts + P_hat;
J=[-l1*sin(q1) -l2*sin(q2); l1*cos(q1) l2*cos(q2)];
Te = Ki*(Pot-P_hat);
Fe=inv(transpose(J))*Te;
end

function [yd,dyd,ddyd]= trajectory_planner(t0,tf,yi,dyi,yf,dyf,z)
global m1 m2 ts Ki kp kd ts l1 l2 l0 
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
yd=a0+(a1*t)+(a2*t*t)+(a3*t*t*t);
dyd=(a1)+(2*a2*t)+(3*a3*t*t);
ddyd=2*a2+ (6*a3*t);
end

function [q,dq]=inverse_kine(xd,dyd,yd)
global m1 m2 ts Ki kp kd ts l1 l2 l0 
% alpha=acos((l1^2+l2^2-(yd^2+z^2))/2*l1*l2);
% theta2=(pi-alpha);
% theta1=atan(z/yd)+atan(l2*sin(theta2)/(l1+l2*cos(theta2)));
% theta2=theta1-theta2;
% q=[theta1; theta2];
syms q1 q2
eq1=l1*cos(q1)+l2*cos(q2)==xd;
eq2=l1*sin(q1)+l2*sin(q2)==yd;
res = solve(eq1,eq2);
q=[double(res.q1(1)); double(res.q2(1))];
theta1=double(res.q1(1));
theta2=double(res.q2(1));
J=[-l1*sin(theta1) -l2*sin(theta2);
    l1*cos(theta1) l2*cos(theta2)];
dq=inv(J)*[dyd; 0];
end
