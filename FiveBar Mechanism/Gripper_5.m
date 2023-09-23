clc 
clear
l2 = 0.5 ;l3 = 0.3 ; l4 = ((2*2.5^2)^0.5)/10 ;l5 = l4 ;off = 0.3;
g = 0;
m2 = 0.5; m3 = 0.3; m4 = 0.25; m5 = 0.25;
k = 0.8;  d = (l2^2+l3^2)^0.5;
q20 = pi/2; q30 =0; q40 = pi/2 + pi/4; q50 = pi/4;
syms q2(t) q3(t) q4(t) q5(t) lam1 lam2 ext
syms dq2 dq3 dq4 dq5 DDq2 DDq3 DDq4 DDq5
Dq2 = diff(q2); Dq3 = diff(q3); Dq4 = diff(q4); Dq5 = diff(q5);
v2sq = (0.5*l2*Dq2)^2; v3sq = (l2*sin(q2)*Dq2 + 0.5*l3*sin(q3)*Dq3)^2 + (l2*cos(q2)*Dq2 + 0.5*l3*cos(q3)*Dq3)^2; v4sq = (l5*sin(q5)*Dq5 + 0.5*l4*sin(q4)*Dq4)^2 + (l5*cos(q5)*Dq5 + 0.5*l4*cos(q4)*Dq4)^2; v5sq = (0.5*l5*Dq5)^2;
i2 = (m2*l2^2)/12 ; i3 = (m3*l3^2)/12 ;i4 = (m4*l4^2)/12 ; i5 = (m5*l5^2)/12 ;
x3 = l2*cos(q2) + 0.5*l3*cos(q3); y3 = l2*sin(q2) + 0.5*l3*sin(q3);

%% Trajectory details
dt=0.01;%0.01
to=0;
tf=2;
timerx = to:dt:tf;
qp2i = pi/2;
qp2di = 0;
qp2f = pi/2 + pi/4;
qp2df= 0;
%% 


%% Trajectory Planning for end-effector 
[theta2,th2dotd,th2ddotd]= trajectory_planner(to,tf,qp2i,qp2di,qp2f,qp2df);
% [yd,ydd,yddd]= trajectory_planner(to,tf,yi,yid,yf,yfd);

for i=1:length(timerx)
    th2(i) = double(subs(theta2,t,timerx(i)));
    th2d(i) = double(subs(th2dotd,t,timerx(i)));
end

disp("Trajectory planned");
%% 

T = 0.5*m2*v2sq + 0.5*i2*Dq2^2 + 0.5*m3*v3sq + 0.5*i3*Dq3^2 + 0.5*m4*v4sq + 0.5*i4*Dq4^2 + 0.5*m5*v5sq + 0.5*i5*Dq5^2;
V = 0.5*m2*g*l2*sin(q2) + m3*g*(l2*sin(q2)+ 0.5*l3*sin(q3)) + m4*g*(l5*sin(q5)+ 0.5*l4*sin(q4)) + 0.5*m5*g*l5*sin(q5) + 0.5*k*((x3^2+y3^2)^0.5 - d)^2;
L = T-V;
cons1 = l2*cos(q2) + l3*cos(q3) - l4*cos(q4) - l5*cos(q5) - off ;
cons2 = l2*sin(q2) + l3*sin(q3) - l4*sin(q4) - l5*sin(q5);

eq1 = diff(diff(L,Dq2)) - diff(L,q2) == ext + lam1* diff(cons1,q2) + lam2* diff(cons2,q2);% add input torque
eq2 = diff(diff(L,Dq3)) - diff(L,q3) == lam1* diff(cons1,q3) + lam2* diff(cons2,q3);
eq3 = diff(diff(L,Dq4)) - diff(L,q4) == lam1* diff(cons1,q4) + lam2* diff(cons2,q4);
eq4 = diff(diff(L,Dq5)) - diff(L,q5) == lam1* diff(cons1,q5) + lam2* diff(cons2,q5);
eq5 = diff(diff(cons1));
eq6 = diff(diff(cons2));
E11 = diff(L,Dq2)*Dq2; E22 = diff(L,Dq3)*Dq3; E33 = diff(L,Dq4)*Dq4; E44 = diff(L,Dq5)*Dq5; 
E=E11+E22+E33+E44-L;

eq1= subs(eq1,[diff(q2(t), t, t) diff(q3(t), t, t) diff(q4(t), t, t) diff(q5(t), t, t)],[DDq2 DDq3 DDq4 DDq5]);
eq2= subs(eq2,[diff(q2(t), t, t) diff(q3(t), t, t) diff(q4(t), t, t) diff(q5(t), t, t)],[DDq2 DDq3 DDq4 DDq5]);
eq3= subs(eq3,[diff(q2(t), t, t) diff(q3(t), t, t) diff(q4(t), t, t) diff(q5(t), t, t)],[DDq2 DDq3 DDq4 DDq5]);
eq4= subs(eq4,[diff(q2(t), t, t) diff(q3(t), t, t) diff(q4(t), t, t) diff(q5(t), t, t)],[DDq2 DDq3 DDq4 DDq5]);
eq5= subs(eq5,[diff(q2(t), t, t) diff(q3(t), t, t) diff(q4(t), t, t) diff(q5(t), t, t)],[DDq2 DDq3 DDq4 DDq5]);
eq6= subs(eq6,[diff(q2(t), t, t) diff(q3(t), t, t) diff(q4(t), t, t) diff(q5(t), t, t)],[DDq2 DDq3 DDq4 DDq5]);
E= subs(E,[diff(q2(t), t), diff(q3(t), t), diff(q4(t), t), diff(q5(t), t)],[dq2 dq3 dq4 dq5]);

[A,B]=equationsToMatrix([eq1, eq2, eq3, eq4, eq5, eq6], [DDq2 DDq3 DDq4 DDq5 lam1 lam2]);
% Dq2 = diff(q2); Dq3 = diff(q3); Dq4 = diff(q4); Dq5 = diff(q5);
% DDq2 = diff(Dq2); DDq3 = diff(Dq3); DDq4 = diff(Dq4); DDq5 = diff(Dq5);
B = subs(B,[diff(q2(t), t), diff(q3(t), t), diff(q4(t), t), diff(q5(t), t)],[dq2 dq3 dq4 dq5]);
q = [q20;q30;q40;q50;0;0] ; dotq = [0;0;0;0;0;0];
tau = [0;]; arrE=[];
disp(B)
kp = 10; kd = 1;
for i=1:length(timerx)
    taui = - kp*(q(1,i)-th2(i)) - kd*(dotq(1,i)-th2d(i))
    AA = subs(A,[q2(t),q3(t),q4(t),q5(t),ext],[q(1,i) q(2,i) q(3,i) q(4,i) taui]);
%     BB = subs(B,[q2(t),q3(t),q4(t), q5(t),  diff(q2(t), t), diff(q3(t), t), diff(q4(t), t), diff(q5(t), t)],[q(1,i) q(2,i) q(3,i) q(4,i) dotq(1,i) dotq(2,i) dotq(3,i) dotq(4,i)]);
    BB = subs(B,[q2(t),q3(t),q4(t), q5(t), dq2, dq3, dq4, dq5, ext ],[q(1,i) q(2,i) q(3,i) q(4,i) dotq(1,i) dotq(2,i) dotq(3,i) dotq(4,i) taui]);
    EE = subs(E,[q2(t),q3(t),q4(t), q5(t), dq2, dq3, dq4, dq5],[q(1,i) q(2,i) q(3,i) q(4,i) dotq(1,i) dotq(2,i) dotq(3,i) dotq(4,i)]);
    ddotq = inv(AA)*BB;

    dotq(:,i+1) = dotq(:,i) + ddotq*dt ;
    q(:,i+1) = q(:,i) + dotq(:,i+1)*dt ;
    tau(:,i+1) = tau(:,i) + taui ;
    arrE(end+1)=EE;
    i
    % q(1,i)
    % th2(i)
end
%%
%Animating
% fig1 = figure() ;
% for i = 1:1:length(timerx)%length(timerx)
%     fk1 = [l2*cos(q(1,i));l2*sin(q(1,i))] ; 
%     fk2 = [l2*cos(q(1,i))+l3*cos(q(2,i));l2*sin(q(1,i))+l3*sin(q(2,i))] ;
%     fk3 = [l5*cos(q(4,i))+l4*cos(q(3,i)) + off;l5*sin(q(4,i))+l4*sin(q(3,i))] ;
%     fk4 = [l5*cos(q(4,i))+ off;l5*sin(q(4,i))]
% 
%     clf ;
%     plot([0,fk1(1,:)],[0,fk1(2,:)],'k','Linewidth',3) ;
%     hold on ;
%     plot([fk1(1,:),fk2(1,:)],[fk1(2,:), fk2(2,:)],'k','Linewidth',3) ;
%     hold on ;
% %     plot([fk2(1,:), fk3(1,:)],[fk2(2,:), fk3(2,:)],'k','Linewidth',3) ;
% %     hold on ;
%     plot([fk2(1,:), fk4(1,:)],[fk2(2,:), fk4(2,:)],'k','Linewidth',3) ;
%     hold on ;
%     plot([off,fk4(1,:)],[0,fk4(2,:)],'k','Linewidth',3) ;
%     hold on ;
%     axis square;
%     xlim([-2 2])
%     ylim([-2 2])
%     drawnow;
% end
fig2 = figure() ;
plot([0,timerx],q(1,:),'r','Linewidth',3) ;
hold on;
plot(timerx,th2,'k','Linewidth',3) ;

fig3 = figure() ;
grid on;
plot(timerx,tau(1,2:202),'k','Linewidth',3) ;

fig4 = figure();
plot(arrE);

function [theta_d,thetadot_d,thetaddot_d]= trajectory_planner(to,tf,thetai,thetadi,thetaf,thetadf)
syms t
Q = [thetai;thetadi;thetaf;thetadf];
t0=to;
B = [1 t0 t0*t0 t0*t0*t0;
    0 1 2*t0 3*t0*t0;
    1 tf tf*tf tf*tf*tf;
    0 1 2*tf 3*tf*tf];
Binv1= inv(B);
A1= Binv1*Q;
a0=A1(1,1);
a1=A1(2,1);
a2=A1(3,1);
a3=A1(4,1);
theta_d=a0+(a1*t)+(a2*t*t)+(a3*t*t*t); % Thita(t)
thetadot_d=(a1)+(2*a2*t)+(3*a3*t*t);
thetaddot_d=2*a2+ (6*a3*t);
end
