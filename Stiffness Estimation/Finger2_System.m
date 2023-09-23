clear; clc;
close all;

global m1 m2 m3 m4 ts Ki kp kd l1 l2 l3 l4 l0 lc1 lc2 lc3 lc4 I1 I2 I3 I4 g stiff ki lambda
ts=0.01; t0=0; tf=2;
tseries = t0:ts:tf;
m1=2; m2=2; m3=2; m4=2; l0=20; Ki=[25 0 0 0; 0 25 0 0; 0 0 25 0; 0 0 0 25]; kp=450; kd=10; ki=0; stiff=5; lambda =0.1;

%%Finger Dimension
l0=2; l1=1.5; l2=1.5; l3=1.5; l4=1.5; lc1=l1/2; lc2=l2/2; lc3=l3/2; lc4=l4/2; I1=l1^2/3; I2=l2^2/3; I3=l3^2/3; I4=l4^2/3; g=9.81;

syms t
x1min=0; x1max=0.5; y1f=2; x2min=l0; x2max=l0-0.5; 
[xd,dxd,ddxd]= trajectory_planner(t0,tf,x1min,0,x1max,0,y1f);
[x2d,dx2d,ddx2d]= trajectory_planner(t0,tf,x2min,0,x2max,0,y1f);
q=[]; dq=[]; prevq=[0;0;0;0]; ddqd=[];
x1f=[]; dx1f=[]; x2f=[]; dx2f=[];
%[y2d,dy2d,ddy2d]= trajectory_planner(t0,tf,y2i,0,y2f,0,z);
for i=1:1:length(tseries)
    i
    x1f(i) = double(subs(xd,t,tseries(i)));
    dx1f(i) = double(subs(dxd,t,tseries(i)));
    x2f(i) = double(subs(x2d,t,tseries(i)));
    dx2f(i) = double(subs(dx2d,t,tseries(i)));
    [x,dx]=inverse_kine(x1f(i),dx1f(i),y1f,x2f(i),dx2f(i),y1f);
    ddqd(:,end+1)=[dx(1)-prevq(1); dx(2)-prevq(2); dx(3)-prevq(3); dx(4)-prevq(4)]./ts;
    q(:,end+1)=[x(1); x(2); x(3); x(4)];
    dq(:,end+1)=[dx(1); dx(2); dx(3); dx(4)];
    prevq(:,end+1)=[dx(1); dx(2); dx(3); dx(4)];
end


T=[0;0;0;0]; Fd=[zeros(1,length(tseries));zeros(1,length(tseries));zeros(1,length(tseries));zeros(1,length(tseries))]; 
Te=[0;0;0;0]; y=q(:,1); dy=dq(:,1);%Fext=[stiff*y1f;zeros(1,length(tseries))];
P_hat=[0;0;0;0]; estF=[]; acty=[]; err=[0;0;0;0]; torq=[];
stiff_est=[0]'; xf=[]; yf=[]; x1i=x1min; y1i=0; x2i=x2min; y2i=0;
Pest =inv(lambda*eye(1)); Kf=[];
for i=1:length(tseries)
 [T,y,dy,Text,err]=system(T,Fd(i),Te,stiff*[-(x1i+l0-x2i);y1i-y1f;x1i+l0-x2i;y2i-y1f],q(:,i),dq(:,i),ddqd(:,i),y,dy,err);
 [Fe,Te,P_hat]=dob(P_hat,T,y,dy);
 [x1i,y1i,x2i,y2i]=forward_kin(y);
 [y_est,stiff_est,e,Pest]=estimateK(-(x1i+l0-x2i),Fe(1),stiff_est,Pest);
 estF(:,end+1)=Fe;
 torq(:,end+1)=stiff*[-(x1i+l0-x2i),y1i-y1f,x1i+l0-x2i,y2i-y1f];%Fext(:,i);
 acty(:,end+1)=y;
 Kf(:,end+1)=stiff_est;
 xf(:,end+1)=[x1i,x2i];
 yf(:,end+1)=[y1i,y2i];
end


fig1=figure();
subplot(4,1,1);plot(q(1,:));hold on;plot(acty(1,:));
legend('Desired Q1','Actual Q1');
title("Theta1 Trajectory Tracking");
xlabel("Time(s)");
ylabel("Angle(rad)");
subplot(4,1,2);plot(q(2,:));hold on;plot(acty(2,:));
legend('Desired Q2','Actual Q2');
title("Theta2 Trajectory Tracking");
xlabel("Time(s)");
ylabel("Angle(rad)");
subplot(4,1,3);plot(q(3,:));hold on;plot(acty(3,:));
legend('Desired Q3','Actual Q3');
title("Theta3 Trajectory Tracking");
xlabel("Time(s)");
ylabel("Angle(rad)");
subplot(4,1,4);plot(q(4,:));hold on;plot(acty(4,:));
legend('Desired Q4','Actual Q4');
title("Theta4 Trajectory Tracking");
xlabel("Time(s)");
ylabel("Angle(rad)");

fig4=figure();
subplot(4,1,1); plot(xf(1,:)); hold on; plot(x1f);
legend('Actual Displacement','Desired Displacement');
title("Trajectory Tracking along x-axis (Mass 1)");
xlabel("Time(s)");
ylabel("Displacement(mm)");
subplot(4,1,2); plot(yf(1,:)); hold on; plot(y1f+zeros(length(tseries))); %
legend('Actual Displacement','Desired Displacement');
title("Trajectory Tracking along y-axis (Mass 1)");
xlabel("Time(s)");
ylabel("Displacement(mm)");
subplot(4,1,3); plot(xf(2,:)); hold on; plot(x2f);
legend('Actual Displacement','Desired Displacement');
title("Trajectory Tracking along x-axis (Mass 2)");
xlabel("Time(s)");
ylabel("Displacement(mm)");
subplot(4,1,4); plot(yf(2,:)); hold on; plot(y1f+zeros(length(tseries))); %
legend('Actual Displacement','Desired Displacement');
title("Trajectory Tracking along y-axis (Mass 2)");
xlabel("Time(s)");
ylabel("Displacement(mm)");

fig2=figure();
subplot(4,1,1);plot(torq(1,:));hold on;plot(estF(1,:));
legend('Actual Force','Estimated Force');
title("Force Estimation along x-axis (Mass 1)");
xlabel("Time(s)");
ylabel("Force(N)");
subplot(4,1,2);plot(torq(2,:));hold on;plot(estF(2,:));
legend('Actual Force','Estimated Force');
title("Force Estimation along y-axis (Mass 1)");
xlabel("Time(s)");
ylabel("Force(N)");
subplot(4,1,3);plot(torq(3,:));hold on;plot(estF(3,:));
legend('Actual Force','Estimated Force');
title("Force Estimation along x-axis (Mass 2)");
xlabel("Time(s)");
ylabel("Force(N)");
subplot(4,1,4);plot(torq(4,:));hold on;plot(estF(4,:));
legend('Actual Force','Estimated Force');
title("Force Estimation along y-axis (Mass 2)");
xlabel("Time(s)");
ylabel("Force(N)");

fig3=figure();
plot(xf(1,:)+l0-xf(2,:),Kf(1,:)); hold on; plot(xf(1,:)+l0-xf(2,:),stiff+zeros(length(tseries)));
legend('Estimated Stiffness','Actual Stiffness');
title("Stiffness Estimation");
xlabel("Displacement(mm)");
ylabel("Stiffness(N/mm)");

function [x1,z1,x2,z2]=forward_kin(y)
global m1 m2 m3 m4 ts Ki kp kd l1 l2 l3 l4 l0 lc1 lc2 lc3 lc4 I1 I2 I3 I4 g stiff ki lambda
q1i=y(1); q2i=y(2); q3i=y(3); q4i=y(4); 
x1 = l1*cos(q1i) + l2*cos(q2i);
z1 = +l1*sin(q1i) + l2*sin(q2i);
x2 = l0+l3*cos(q3i) + l4*cos(q4i);
z2 = l3*sin(q3i) + l4*sin(q4i);
end

function [q,dq]=inverse_kine(x1d,dy1d,y1d,x2d,dy2d,y2d)
global m1 m2 m3 m4 ts Ki kp kd l1 l2 l3 l4 l0 lc1 lc2 lc3 lc4 I1 I2 I3 I4 g stiff ki lambda

syms q1 q2 q3 q4
eq1=l1*cos(q1)+l2*cos(q2)==x1d;
eq2=l1*sin(q1)+l2*sin(q2)==y1d;
res = solve(eq1,eq2);
theta1=double(res.q1(1));
theta2=double(res.q2(1));
eq3=l3*cos(q3)+l4*cos(q4)==x2d-l0;
eq4=l1*sin(q3)+l2*sin(q4)==y2d;
res2 = solve(eq3,eq4);

q=[double(res.q1(1)); double(res.q2(1)); double(res2.q3(2)); double(res2.q4(2))];
J=[-l1*sin(q(1)) -l2*sin(q(2)) 0 0; l1*cos(q(1)) l2*cos(q(2)) 0 0; 0 0 -l3*sin(q(3)) -l4*sin(q(4)); 0 0 l3*cos(q(3)) l4*cos(q(4))];
dq=inv(J)*[dy1d; 0; dy2d; 0];
end

function [y_est,stiffi,e,P]=estimateK(x,y,stiffi,P)
global m1 m2 m3 m4 ts Ki kp kd l1 l2 l3 l4 l0 lc1 lc2 lc3 lc4 I1 I2 I3 I4 g stiff ki lambda
X = [x]';
K = (P*X) / (lambda + X'*P*X);
% K = (P*X)*inv((lambda + X'*P*X));
y_est = X' * stiffi;
e = y-y_est;
stiffi =stiffi + K*e;
P = P - K * X' * P;
end

function [T,y,dy,Text,err]=system(T,Fd,Te,Fext,yd,dyd,ddqd,y,dy,err)
global m1 m2 m3 m4 ts Ki kp kd l1 l2 l3 l4 l0 lc1 lc2 lc3 lc4 I1 I2 I3 I4 g stiff ki lambda
q1=y(1); q2=y(2); q3=y(3); q4=y(4); 

%% Dynamics
M = [(m1*l1^2/4)+m2*l1^2+I1, (m2*l1*l2*cos(q1-q2))/2, 0, 0;
          (m2*l1*l2*cos(q1-q2))/2, (m2*l2^2)/4 + I2, 0, 0;
          0, 0, (m3*l3^2/4)+m4*l3^2+I3, (m4*l3*l4*cos(q3-q4))/2;
          0, 0, (m4*l3*l4*cos(q3-q4))/2, (m4*l4^2)/4 + I4] ;  
C = [0,(m2*l1*l2*sin(q2-q1))/2, 0, 0;
     -(m2*l1*l2*sin(q2-q1))/2, 0, 0, 0;
     0, 0, 0,(m4*l3*l4*sin(q4-q3))/2;
     0, 0, -(m4*l3*l4*sin(q4-q3))/2, 0] ;
G=0; 
J=[-l1*sin(q1) -l2*sin(q2) 0 0; l1*cos(q1) l2*cos(q2) 0 0; 0 0 -l3*sin(q3) -l4*sin(q4); 0 0 l3*cos(q3) l4*cos(q4)];
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
global m1 m2 m3 m4 ts Ki kp kd l1 l2 l3 l4 l0 lc1 lc2 lc3 lc4 I1 I2 I3 I4 g stiff ki lambda
q1=y(1); q2=y(2); q3=y(3); q4=y(4); 
M = [(m1*l1^2/4)+m2*l1^2+I1, (m2*l1*l2*cos(q1-q2))/2, 0, 0;
          (m2*l1*l2*cos(q1-q2))/2, (m2*l2^2)/4 + I2, 0, 0;
          0, 0, (m3*l3^2/4)+m4*l3^2+I3, (m4*l3*l4*cos(q3-q4))/2;
          0, 0, (m4*l3*l4*cos(q3-q4))/2, (m4*l4^2)/4 + I4] ;  
C = [0,(m2*l1*l2*sin(q2-q1))/2, 0, 0;
     -(m2*l1*l2*sin(q2-q1))/2, 0, 0, 0;
     0, 0, 0,(m4*l3*l4*sin(q4-q3))/2;
     0, 0, -(m4*l3*l4*sin(q4-q3))/2, 0] ;
Pot = M*dy;
err = Pot-P_hat;
P_hat= (T+C'*dy+Ki*err)*ts + P_hat;
J=[-l1*sin(q1) -l2*sin(q2) 0 0; l1*cos(q1) l2*cos(q2) 0 0; 0 0 -l3*sin(q3) -l4*sin(q4); 0 0 l3*cos(q3) l4*cos(q4)];
Te = Ki*(Pot-P_hat);
Fe=inv(transpose(J))*Te;
end

function [yd,dyd,ddyd]= trajectory_planner(t0,tf,yi,dyi,yf,dyf,z)
global m1 m2 m3 m4 ts Ki kp kd l1 l2 l3 l4 l0 lc1 lc2 lc3 lc4 I1 I2 I3 I4 g stiff ki lambda
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

