clear; clc;

m1=0.4;m2=0.4;l1=0.5;l2=0.5;l1g=0.5*l1;l2g=0.5*l2;j1=m1*l1^2/12;j2=m2*l2^2/12;g=9.81;
ts=0.01;
tseries = 0:0.01:1;

yd = 0:0.01:1;
%1<1-CBr
r=1/29889;%0.739;%0.00055;
I=[1 0; 0 1]; dets=[];

xi=[pi/2; 0];
x0=[pi/2; 0]; dx0=[0; 0]; ddx0=[0; 0]; error0=[(yd(1)-x0(1)); (yd(2)-x0(2))];u0=[0; 0];
q1=[];q2=[];
x1=x0(1); x2=x0(2); dx=dx0; error=error0; derror=error0;ierror=error0;

uprev=[];q1s=[];
for i=1:100
%%System Dynamics
x1=x0(1); x2=x0(2);
% M=[m1*l1g^2+m2*l1^2+j1 m2*l1*l2g*cos(x1-x2);
%     m2*l1*l2g*cos(x1-x2) m2*l2g^2+j2];
% C=-m2*l1*l2g*g*sin(x1-x2);
% G=[(m1*l1g+m2*l1)*g*cos(x1); m2*l2g*g*cos(x2)];


%%Controller
% a=0.01; k=65; r=k*inv(M);
% u = a*(u0+r*(error));

%% D-Type
% u=(u0)+(r*inv(M))*abs(derror);
% dets(end+1)=det(I-inv(M)*r);%to check convergence

%% P-Type
% u=u0+(r*inv(M))*abs(error);

%% PID Type
kp=1;kd=1;ki=1;
u=u0+kp*(error)+ki*ierror+kd*derror;

uprev(:,end+1)=u;
%%Solving the system
ddx = inv(M)*(u-C*dx0-G);%ddx0 + 
dx = dx0+ddx*ts;
x = x0+dx*ts;
ddx0=ddx; dx0=dx; x0=x;  derror = [((yd(i)-x(1))-error(1))/ts; ((yd(i)-x(2))-error(2))/ts]; ierror= error+[(yd(i)-x(1)); (yd(i)-x(2))];error=[(yd(i)-x(1)); (yd(i)-x(2))]; u0=u; 
q1(end+1)=x(1);q2(end+1)=x(2);
end
q1s(end+1,:)=q1; q1=[];

for j=1:5
    u0=uprev;
    uprev=[];
    x0=[pi/2; 0];
for i=1:100
%%System Dynamics
x1=x0(1); x2=x0(2);
M=[m1*l1g^2+m2*l1^2+j1 m2*l1*l2g*cos(x1-x2);
    m2*l1*l2g*cos(x1-x2) m2*l2g^2+j2];
C=-m2*l1*l2g*g*sin(x1-x2);
G=[(m1*l1g+m2*l1)*g*cos(x1); m2*l2g*g*cos(x2)];

%%Controller
% a=0.01; k=65; r=k*inv(M);
% u = a*(u0+r*(error));
%% D-Type
% u=(u0(i))+(r*inv(M))*abs(derror);
% dets(end+1)=abs(det(I-inv(M)*r));%to check convergence

%% P-Type
% u=u0(i)+(r*inv(M))*abs(error);

%% PID Type
kp=1;kd=1;ki=1;
u=u0(i)+kp*(error)+ki*ierror+kd*derror;

uprev(:,end+1)=u;
%%Solving the system
ddx = inv(M)*(u-C*dx0-G);%ddx0 + 
dx = dx0+ddx*ts;
x = x0+dx*ts;
ddx0=ddx; dx0=dx; x0=x;  derror = [((yd(i)-x(1))-error(1))/ts; ((yd(i)-x(2))-error(2))/ts]; ierror= [(yd(i)-x(1))+error(1); (yd(i)-x(2))+error(2)]; error=[(yd(i)-x(1)); (yd(i)-x(2))]; %u0=u; 
q1(end+1)=x(1);q2(end+1)=x(2);
end
q1s(end+1,:)=q1;
q1=[];
end

disp(max(dets));
disp(min(dets));
plot(yd);
hold on;
plot(q1s(end,:));
hold off;