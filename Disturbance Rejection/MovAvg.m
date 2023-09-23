clc; clear all; close all;

t0=0; ts=0.1; tf=5;
t=t0:ts:tf;
a0=5; w=0.5;

%System Parameters 
a1=1; a2=0.7; b0=1; b1=0.7; c1=1; c2=0.5;
r1=((a2-c2)*b0*b1+(c1-a1)*b1^2)/(b1^2+a1*b0*b1+a2*b0^2);
s0=(b1*(a1^2-a2-c1*a1+c2)+(b0*c1*a2-b0*a1*a2))/(b1^2+a1*b0*b1+a2*b0^2);
s1=(b1*(a1*a2-c1*a2)+(b0*a2*c2-b0*a2^2))/(b1^2+a1*b0*b1+a2*b0^2);
u_t=zeros(length(t));% a0*sin(w*t);%
e=rand(length(t),1)';

y_t=[0 0]; uc_t=[0 0]; err=[0 0];
for i=3:length(t)
y_t(end+1)=b0*uc_t(i-1)+b1*u_t(i-2)+e(i)+c1*e(i-1)+c2*e(i-2)-a1*y_t(i-1)-a2*y_t(i-2);
%y_t(end+1)=b0*uc_t(i-1)+b1*uc_t(i-2)+e(i)-a1*y_t(i-1)-a2*y_t(i-2);
uc_t(end+1)=-r1*u_t(i-1)-s0*y_t(i)-s1*y_t(i-1);
err(end+1)=e(i);%+r1*e(i-1);
end

fig1=figure();
plot(u_t);
fig2=figure();
plot(y_t);
hold on; 
plot(err);
legend('Output Signal','Error');