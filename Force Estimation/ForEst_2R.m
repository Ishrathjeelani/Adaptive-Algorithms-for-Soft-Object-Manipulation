clear; clc;
close all;

global m1 m2 l1 l2 l1g l2g j1 j2 g ts Ki Mi
m1=0.8;m2=0.4;l1=0.8;l2=0.4;l1g=0.5*l1;l2g=0.5*l2;j1=m1*l1^2/12;j2=m2*l2^2/12;g=9.81;
ts=0.01;
tseries = 0:0.01:50;


tm=[2+zeros( 1, length(tseries));1+zeros( 1, length(tseries))];%[0.15 + zeros( 1, length(tseries));0.15 + zeros( 1, length(tseries))];
text=[ 5+zeros( 1, length(tseries));1+zeros( 1, length(tseries))];%[0.5 + zeros( 1, length(tseries)); 0.5 + zeros( 1, length(tseries))];
x=[0; pi/4]; dx=[0; 0]; r=tm(:,1); t_dis=[r]; p_cap=0; p=0; Mi=[m1 0;0 m2]; Ki=[0.2 0; 0 0.2];
theta1=[]; theta2=[];
for i=1:length(tseries)
[q,dq]=robot(tm(i),text(i),x,dx);
x=q;dx=dq;
[p_prev,p_cap_prev,r_prev]=dob(x,dx,tm(i),text(i),p,p_cap,r);
p_cap=p_cap_prev; r=r_prev;p=p_prev;
t_dis(:,end+1)=r;
theta1(end+1)=q(1);
theta2(end+1)=q(2);
end

plot(t_dis(2,:));
title("Torque vs Time");
xlabel("Time(s)");
ylabel("Torque(N-m)");
hold on;
plot(t_dis(1,:));
plot(text(1,:));
plot(text(2,:));
legend('Disturbance Joint1','Disturbance Joint2','Applied Torque1','Applied Torque2');

fig2=figure();
plot(theta1);
hold on;
plot(theta2);
legend('Joint1','Joint2');

function [x,dx]=robot(tj,text,x,dx)
global m1 m2 l1 l2 l1g l2g j1 j2 g ts Ki
x1=x(1); x2=x(2);
M=[m1*l1g^2+m2*l1^2+j1 m2*l1*l2g*cos(x1-x2);
    m2*l1*l2g*cos(x1-x2) m2*l2g^2+j2];
C=-m2*l1*l2g*g*sin(x1-x2);
G=[(m1*l1g+m2*l1)*g*cos(x1); m2*l2g*g*cos(x2)];
u=tj+text;
ddx = inv(M)*(u-C*dx-G);
dx = dx+ddx*ts;
x = x+dx*ts;
%x=rem(x,2*pi);
end

function [p,p_cap,r]=dob(x,dq,tj,t_dis,p,p_cap,r)
global m1 m2 l1 l2 l1g l2g j1 j2 g ts Ki Mi
x1=x(1); x2=x(2); %tj=tj+t_dis;
M=[m1*l1g^2+m2*l1^2+j1 m2*l1*l2g*cos(x1-x2);
    m2*l1*l2g*cos(x1-x2) m2*l2g^2+j2];
C=[-m2*l1*l2g*g*sin(x1-x2)];%
G=[(m1*l1g+m2*l1)*g*cos(x1); m2*l2g*g*cos(x2)];
p=p+M*dq*ts; dp_cap=tj+transpose(C).*dq-G+r;
p_cap=p_cap+dp_cap*ts;
%dp=tau+transpose(C)*dq-G+Ki(p-p_cap);
r=Ki*(p-p_cap);
end