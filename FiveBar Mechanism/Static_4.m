% clc 
% clear
l2 = 0.5 ;l3 = 0.3 ; l4 = 0.5 ;;off = 0.3; fin = 0.3;
m2 = 0.5; m3 = 0.3; m4 = 0.25; 
g=0;
i2 = 0.5*m2*l2^2/12; i3 = 0.5*m3*l3^2/12; i4 = 0.5*m4*l4^2/12;
mass = 5; gg= 9.81;
umu = 0.3; w= mass*gg;
Fx = w/umu; Fy = -w; Text = -0.5*fin*Fx + 0.5*l3*Fy;

%%kinematics
syms q2 q3 q4 %q2(t) q3(t) q4(t)
q20 = pi/6;
x2x = l2*cos(q2); x2y = l2*sin(q2); x3x = l3*cos(q3); x3y = l3*sin(q3); x4x = -l4*cos(q4); x4y = -l4*sin(q4); x1x = -off; x1y = 0; 
% A = [x2x x3x x4x x1x;
%     x2y x3y x4y x1y];
%sol = solve(x2x + x3x + x4x + x1x == 0, x2y + x3y + x4y + x1y == 0,q3,q4);
sol = solve(l2*cos(q20)+l3*cos(q3)-l4*cos(q4)-off==0,l2*sin(q20)+l3*sin(q3)-l4*sin(q4)==0,q3,q4);
% sol = solve(x2x+x3x+x4x+x1x==0,x2y+x3y+x4y==0,q3,q4);
val_q3 = double(subs(sol.q3,[q2],[q20]));val_q4 = double(subs(sol.q4,[q2],[q20]));
q30 = val_q3(1); q40 = val_q4(1);
syms w2 w3 w4
w20 = 0.2;
v2y = w20*diff(x2x); v2x = w20*diff(x2y); v3y = w3*diff(x3x); v3x = w3*diff(x3y); v4y = -w4*diff(x4x); v4x = -w4*diff(x4y); 
vel_sol = solve(v2x+v3x+v4x==0,v2y+v3y+v4y==0,w3,w4);
val_w3 = double(subs(vel_sol.w3,[q2,q3,q4,w2],[q20,q30,q40,w20]));val_w4 = double(subs(vel_sol.w4,[q2,q3,q4,w2],[q20,q30,q40,w20]));
w30 = val_w3(1); w40 = val_w4(1);
syms a2 a3 a4
a20 = 0;
acc_sol = solve(-a20*l2*sin(q20)-l2*w20^2*cos(q20)-a3*l3*sin(q30)-l3*w30^2*cos(q30)+a4*l4*sin(q40)+l4*w40^2*cos(q40)==0, a20*l2*cos(q20)-l2*w20^2*sin(q20)+a3*l3*cos(q30)-l3*w30^2*sin(q30)-a4*l4*cos(q40)+l4*w40^2*sin(q40)==0, a3, a4);
a30 = double(acc_sol.a3); a40 = double(acc_sol.a4);
a2x = -l2*a20*sin(q20) - l2*w20^2*cos(q20); a2y = l2*a20*cos(q20) - l2*w20^2*sin(q20); a3x = -l2*a20*sin(q20) - l2*w20^2*cos(q20)-l3*a30*sin(q30) - l3*w30^2*cos(q30); a3y = l2*a20*cos(q20) - l2*w20^2*sin(q20)+l3*a30*cos(q30) - l3*w30^2*sin(q30); a4x = -l4*a40*sin(q40) - l4*w40^2*cos(q40); a4y = l4*a40*cos(q40) - l4*w40^2*sin(q40);
a2x=0.5*a2x; a2y=0.5*a2y; a3x=0.5*a3x; a3y=0.5*a3y; a4x=0.5*a4x; a4y=0.5*a4y; 

r2x = 0.5*l2*cos(q20); r2y = 0.5*l2*sin(q20); r3x = 0.5*l3*cos(q30); r3y = 0.5*l3*sin(q30); r4x = 0.5*l4*cos(q40); r4y = 0.5*l4*sin(q40);
alpha2 = a20;alpha3 = a30;alpha4 = a40;
disp(alpha4);
eqns = [F12x+F32x==m2*a2x,F12y+F32y==m2*a2y, -r2y*F12x+r2x*F12y-r2y*F32x+r2x*F32y==i2*alpha2,
   -F32x+F43x==m3*a3x-Fx,-F32y+F43y==m3*a3y-Fy, r3y*F32x-r3x*F43y-r3y*F32x+r3x*F43y==i3*alpha3+Text,
   -F43x+F14x==m4*a4x,-F43y+F14y==m4*a4y, r4y*F43x-r4x*F43y-r4y*F14x+r4x*F14y+Tout==i4*alpha4];
sol = solve(eqns,[F12x F12y F32x F32y F43x F43y F45x F45y F14x F14y Tout]);
disp(double(sol.Tout));