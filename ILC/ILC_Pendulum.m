clear; clc; close all;

m=1;l=1;g=9.81;
ts=0.001;
tseries = 0:0.001:5;
n=120;
% %% Simple PID
% x0=-pi/2; dx0=0; yd=pi/2;
% x=x0; dx=dx0;
% kp=2;kd=0.5;ki=0.02;
% pos=[x]; vel=[dx];error=0;errord=0;
% 
% for i=1:length(tseries)
% M=m*l^2; C=0; G=m*g*l*sin(x);
% u=kp*(yd-x)+kd*(-dx)+ki*(error);
% ddx = (u-C*dx-G)/M;%ddx0 + 
% dx = dx+ddx*ts;
% x = x+dx*ts;
% error=error+(yd-x);
% errord=errord+(-dx);
% pos(end+1)=x;
% vel(end+1)=dx;
% end
% %%

x0=0; dx0=0; yd=linspace(0,pi/4,length(tseries)+1);%yd=linspace(0,pi/2,length(tseries)+1);
x=x0; dx=dx0;
r=0.12; a=1; %1.5*m*l^2;a=0.1;
pos=[x]; vel=[dx];error=0;derror=0;u0=zeros( 1, length(tseries));err=[];e(1,1)=0;
for j=1:n
    x=x0;dx=dx0;error=0;derror=0;errori=0;Fe=0;
for i=1:length(tseries)
    e(j,i) = yd(i)-x;
    M=m*l^2; C=0; G=m*g*l*sin(x);
    % u=u0(1,i)+r*(error)+a*derror;
    ddx = (u0(j,i)-C*dx-G)/M;
    dx = dx+ddx*ts;
    x = x+dx*ts;
    %x=rem(x,2*pi);
    % error=(yd-x);
    u0(j+1,i) = a*u0(j,i)+r*e(j,i);
    pos(end+1)=x;
    vel(end+1)=dx;
end
end
fig=figure();
subplot(2,3,1);
plot(yd);
title("Position vs Time");
xlabel("Time(s)");
ylabel("Angle(rad)");
hold on;
plot(pos(1:length(tseries)));
hold off;
subplot(2,3,2);
plot(yd);
title("Position vs Time");
xlabel("Time(s)");
ylabel("Angle(rad)");
hold on;
plot(pos(length(tseries)*n/2:length(tseries)*(n/2+1)));
hold off;
subplot(2,3,3);
plot(yd);
title("Position vs Time");
xlabel("Time(s)");
ylabel("Angle(rad)");
hold on;
plot(pos(length(tseries)*(n-1):length(tseries)*(n)));
hold off;
subplot(2,1,2);
plot(vel(length(tseries)*(n-1):length(tseries)*(n)));
title("Velocity vs Time");
xlabel("Time(s)");
ylabel("Anglular Velocity(rad/s)");
fig2=figure();
subplot(3,1,1);
plot(e(1,:));
title("Error trial1 vs Time");
xlabel("Time(s)");
ylabel("Error(rad)");
subplot(3,1,2);
plot(e(n/2,:));
title("Error halfway vs Time");
xlabel("Time(s)");
ylabel("Error(rad)");
subplot(3,1,3);
plot(e(n,:));
title("Error last trial vs Time");
xlabel("Time(s)");
ylabel("Error(rad)");