% clc 
% clear
% l2 = 0.047 ;l3 = 0.018 ; l4 = 0.048 ;l5 = l4 ;off = 0.018; fin = 0.3;
l2 = 0.085 ;l3 = 0.029 ; l4 = 0.071 ;l5 = 0.042 ;off = 0.019; fin = 0.026;res=104;
l4=((l5+l4*cos((180-res)*pi/180))^2 + (l4*sin((180-res)*pi/180))^2)^0.5;
m2 = 0.5; m3 = 0.3; m4 = 0.25; m5 = 0.25;
g=0;
i2 = 0.5*m2*l2^2/12; i3 = 0.5*m3*l3^2/12; i4 = 0.5*m4*l4^2/12; i5 = 0.5*m5*l5^2/12;
mass = 0.5; gg= 9.81;
umu = 0.5; w= mass*gg;
Fx = w/umu; Fy = w; Text = -0.5*fin*Fy + 0.5*l3*Fx;

syms q2 q3 q4 
q20 = pi/4;
x2x = l2*cos(q2); x2y = l2*sin(q2); x3x = l3*cos(q3); x3y = l3*sin(q3); x4x = -l4*cos(q4); x4y = -l4*sin(q4); x1x = -off; x1y = 0; 
sol = solve(l2*cos(q20)+l3*cos(q3)-l4*cos(q4)-off==0,l2*sin(q20)+l3*sin(q3)-l4*sin(q4)==0,q3,q4);
val_q3 = double(subs(sol.q3,[q2],[q20]));val_q4 = double(subs(sol.q4,[q2],[q20]));
q30 = val_q3(1); q40 = val_q4(1);
% r2x = 0.5*l2*cos(q20); r2y = 0.5*l2*sin(q20); r3x = 0.5*l3*cos(q30); r3y = 0.5*l3*sin(q30); r4x = 0.5*l4*cos(q40); r4y = 0.5*l4*sin(q40);
r2x = 0.5*l2*cos(q20); r2y = 0.5*l2*sin(q20); r3x = l2*cos(q20)+0.5*l3*cos(q30); r3y = l2*sin(q20)+0.5*l3*sin(q30); r4x = 0.5*l4*cos(q40); r4y = 0.5*l4*sin(q40);
%r2x = 0; r2y = l2/2; r3x = l3/2; r3y = l2/2; r4x = l3 + l4*cos(pi/4)/2; r4y = l2 - l4*sin(pi/4)/2; r5x = off + l5*cos(pi/4); r5y = l5*sin(pi/4);
a2x = 0; a2y =0; a3x = 0; a3y =0; a4x = 0; a4y =0; a5x = 0; a5y = 0; alpha2 = 0;alpha3 = 0;alpha4 = 0;alpha5 = 0;

syms F12x F12y F32x F32y F43x F43y F45x F45y F14x F14y Tout
% eqns = [F12x+F23x == m2*a2x,F12y+F23y == m2*a2y+ m2*g,r2x*F12y-r2y*F12x+r2x*F23y-r2y*F23x == i2*alpha2,
%     -F23x+F34x == -Fx+m3*a3x,-F23y+F34y  == m3*g+Fy+m3*a3y,Text-r3x*F23y+r3y*F23x+r3x*F34y-r3y*F34x == i3*alpha3,
%     -F34x+F45x == m4*a4x,-F34y+F45y  == m4*g+m4*a4y,-r4x*F34y+r4y*F34x+r4x*F45y-r4y*F45x == i4*alpha4,
%     -F45x+F51x == m5*a5x,-F45y+F51y ==
%     m5*g+m5*a5y,-r5x*F45y+r5y*F45x+r5x*F51y-r5y*F51x-Tout==
%     i5*alpha5];%Tinp+ %5bar


% disp(double(A));
% disp(double(B));
%sol = double(A\B)
m=[0.05,0.1,0.15,0.2,0.25];
h=[0.02,0.03,0.04,0.05,0.06];
tau =[0;];
%%working code
% for i=1:5
%     mass = m(i); gg= 9.81;
%     umu = 0.3; w= mass*gg;
%     Fx = w/umu; Fy = w; Text = -0.5*fin*Fx + 0.5*l3*Fy;
%     eqns = [F12x+F32x==m2*a2x,F12y+F32y==m2*a2y, -r2y*F12x+r2x*F12y-r2y*F32x+r2x*F32y+Tout==i2*alpha2,
%    -F32x+F43x==m3*a3x+Fx,-F32y+F43y==m3*a3y-Fy, r3y*F32x-r3x*F43y-r3y*F32x+r3x*F43y==-Text+i3*alpha3,
%    -F43x+F14x==m4*a4x,-F43y+F14y==m4*a4y, r4y*F43x-r4x*F43y-r4y*F14x+r4x*F14y+Tout*(l4*sin(q40)/(l2*sin(q20)))==i4*alpha4];
%     [A,B] = equationsToMatrix(eqns);
%     sol = solve(eqns,[F12x F12y F32x F32y F43x F43y F45x F45y F14x F14y Tout]);
%     taui =double(sol.Tout);
%     tau(:,i+1) = tau(:,i) + taui ;
% end

% Tip =0.02;
% for i=1:5
%     mass = m(i); gg= 9.81;
%     umu = 0.3; w= mass*gg;
%     Fx = w/umu; Fy = w; Text = -0.5*fin*Fx+ 0.5*l3*Fy;
%     hi = h(i);
% 
%     q20 = asin(hi/l2);
%     x2x = l2*cos(q2); x2y = l2*sin(q2); x3x = l3*cos(q3); x3y = l3*sin(q3); x4x = -l4*cos(q4); x4y = -l4*sin(q4); x1x = -off; x1y = 0; 
%     sol = solve(l2*cos(q20)+l3*cos(q3)-l4*cos(q4)-off==0,l2*sin(q20)+l3*sin(q3)-l4*sin(q4)==0,q3,q4);
%     val_q3 = double(subs(sol.q3,[q2],[q20]));val_q4 = double(subs(sol.q4,[q2],[q20]));
%     q30 = val_q3(1); q40 = val_q4(1);
%     r2x = 0.5*l2*cos(q20); r2y = 0.5*l2*sin(q20); r3x = 0.5*l3*cos(q30); r3y = 0.5*l3*sin(q30); r4x = 0.5*l4*cos(q40); r4y = 0.5*l4*sin(q40);
%     syms w2 w3 w4
%     w20 = 0;
%     v2y = w20*diff(x2x); v2x = w20*diff(x2y); v3y = w3*diff(x3x); v3x = w3*diff(x3y); v4y = -w4*diff(x4x); v4x = -w4*diff(x4y); 
%     vel_sol = solve(v2x+v3x+v4x==0,v2y+v3y+v4y==0,w3,w4);
%     val_w3 = double(subs(vel_sol.w3,[q2,q3,q4,w2],[q20,q30,q40,w20]));val_w4 = double(subs(vel_sol.w4,[q2,q3,q4,w2],[q20,q30,q40,w20]));
%     w30 = val_w3(1);
%     w40 = val_w4(1);
%     syms a2 a3 a4
%     a20 = 0;
%     acc_sol = solve(-a20*l2*sin(q20)-l2*w20^2*cos(q20)-a3*l3*sin(q30)-l3*w30^2*cos(q30)+a4*l4*sin(q40)+l4*w40^2*cos(q40)==0, a20*l2*cos(q20)-l2*w20^2*sin(q20)+a3*l3*cos(q30)-l3*w30^2*sin(q30)-a4*l4*cos(q40)+l4*w40^2*sin(q40)==0, a3, a4);
%     a30 = double(acc_sol.a3);
%     a40 = double(acc_sol.a4);
%     alpha2 = a20;alpha3 = a30;alpha4 = a40;a2x = -l2*a20*sin(q20) - l2*w20^2*cos(q20); a2y = l2*a20*cos(q20) - l2*w20^2*sin(q20); a3x = -l2*a20*sin(q20) - l2*w20^2*cos(q20)-l3*a30*sin(q30) - l3*w30^2*cos(q30); a3y = l2*a20*cos(q20) - l2*w20^2*sin(q20)+l3*a30*cos(q30) - l3*w30^2*sin(q30); a4x = -l4*a40*sin(q40) - l4*w40^2*cos(q40); a4y = l4*a40*cos(q40) - l4*w40^2*sin(q40);
%     a2x=0.5*a2x; a2y=0.5*a2y; a3x=0.5*a3x; a3y=0.5*a3y; a4x=0.5*a4x; a4y=0.5*a4y; 
% 
%     eqns = [F12x+F32x+Fx==m2*a2x,F12y+F32y+Fy==m2*a2y, -r2y*F12x+r2x*F12y-r2y*F32x+r2x*F32y==i2*alpha2,
%    -F32x+F43x==m3*a3x,-F32y+F43y==m3*a3y, r3y*F32x-r3x*F43y-r3y*F32x+r3x*F43y==i3*alpha3,
%    -F43x+F14x==m4*a4x,-F43y+F14y==m4*a4y, r4y*F43x-r4x*F43y-r4y*F14x+r4x*F14y+Tout==i4*alpha4];
%     % [A,B] = equationsToMatrix(eqns);
%     % inv(A)
%     sol = solve(eqns,[F12x F12y F32x F32y F43x F43y F45x F45y F14x F14y Tout]);
%     taui =sol.Tout
%     tau(:,i+1) = tau(:,i) + taui ;
% end


 disp(tau);
 fig3 = figure() ;
plot(m,tau(1,2:length(tau)),'k','Linewidth',3) ;