% clear; clc;
% close all;

global Ct Kt Gr Ra Ke a_freq g_freq m l g Ji ts

Kt=1.5; Gr=1.6; Ra=25; Ke=0.15; Ct=0.5;
a_freq=10;%800*10^(3); 
g_freq=10;%900*10^(3);

m=0.4;l=0.5;g=9.81;q=[0;0];dq=[0;0];
ts=0.01;
tseries = 0:0.01:100;
Jacob=[l*cos(q); l*sin(q)]; Ji=0.01; ext=0.25; tm=[zeros( 1, length(tseries)); zeros( 1, length(tseries))];%0.15 + zeros( 1, length(tseries));

text=ext + zeros( 1, length(tseries));
id=tm.*Ra/(Ct*Kt*Gr);
t_dis=[]; t_m=[];

for i=1:length(tseries)-1
ia(:,i)=motor(id(:,i));
tj(:,i)=ia(:,i)*Kt*Gr;
[q,dq]=robot(tj(:,i),text(i),q,dq);
tl(:,i)=dob(ia(:,i),q,dq);
i_dis(:,i)=feedback(tl(:,i));
id(:,i+1)=id(:,i+1)-i_dis(:,i);
t_dis(:,end+1)=tl(:,i);
t_m(:,end+1)=tj(:,i);
end

% plot(text);
title("Disturbance Torque");
xlabel("Time(s)");
ylabel("Torque(N-m)");
hold on;
plot(t_dis(1,:));
plot(t_dis(2,:));
% axis([0 length(tseries) 0 2*ext]);
hold off;

%plot(id);

function ia=motor(id)
global Ct Kt Gr Ra Ke a_freq g_freq m l g Ji ts
ia= id*Ct/Ra;
end

function [x,dx]=robot(tj,text,x,dx)
global Ct Kt Gr Ra Ke a_freq g_freq m l g Ji ts
M=m*l^2; C=0; G=m*g*l*sin(x);
u=tj+text;
ddx = (u-C*dx-G)/M;
dx = dx+ddx*ts;
x = x+dx*ts;
x=rem(x,2*pi);
end

function tl=dob(ia,q,dq)
global Ct Kt Gr Ra Ke a_freq g_freq m l g Ji ts
J=[50 0; 0 50];%[l*cos(q(1)) 0; 0 l*sin(q(2))];
dqf=dq;%*g_freq*exp(-g_freq*ts);
tl=ia*Kt*Gr-J*dqf;%a_freq*exp(-a_freq*ts)*(dqf*a_freq*Ji+ia*Kt*Gr)-a_freq*Ji*dqf;
end

function i_dis=feedback(tl)
global Ct Kt Gr Ra Ke a_freq g_freq m l g Ji ts
i_dis=tl/(Kt*Gr);
end
