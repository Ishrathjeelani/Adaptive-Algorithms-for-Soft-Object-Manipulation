clc; clear all; close all;

% ap = [1 -1.5 0.7];
% bp = [0 1 0.5];
% cp = [ 1 0 0];
% 
% cnum = [1.5 -0.7];
% cden = [1 0.5];
% 
% nytet = conv(cp,cden);
% dytet = conv(ap,cden) + conv(bp,cnum);
% [nytet, dytet] = minreal(nytet,dytet);
% 
% e = normrnd(0,1,100,1);
% 
% yclp = filter(nytet,dytet,e);
% yolp = filter(cp,ap,e);
% 
% plot([yclp yolp]);


t0=0; ts=0.1; tf=5;
t=t0:ts:tf;
a0=2; w=2;
u_t=a0*sin(w*t);%zeros(length(t));% 
e=rand(length(t),1)';

%System Parameters 
a1=1; a2=0.7; b0=1; b1=0.7; c1=0; c2=0;

y_t=[0 0]; uc_t=[0 0]; err=[e(1) e(2)];
for i=3:length(t)
y_t(end+1)=b0*uc_t(i-1)+b1*u_t(i-2)+e(i)+c1*e(i-1)+c2*e(i-2)+a1*y_t(i-1)-a2*y_t(i-2);
uc_t(end+1)=-b1*u_t(i-1)+(c1-a1)*y_t(i)-(c2-a2)*y_t(i-1);
err(end+1)=e(i);
end

fig1=figure();
plot(u_t);
fig2=figure();
plot(y_t);
hold on; 
plot(err);
legend('Output Signal','Error');