% Kellen C. Schroeter
% University of Colorado 2016
% Course Project

%---------------------
%Set Up
tic;
clc;
clear;
close all;

syms t;
time = linspace(0,2*pi,100).';
time(1,:) = .0001;

%---Read the Text File---
Pos_read = table2array(readtable('mav_data_txt.txt'));
%convert data to symbolic function
timeread = Pos_read(:,1);
tread = Pos_read(:,2);
nread = Pos_read(:,3);
bread = Pos_read(:,1);
cvt = coeffvalues(fit(timeread,tread,'poly4'));
cvn = coeffvalues(fit(timeread,nread,'poly7'));
cvb = coeffvalues(fit(timeread,bread,'poly2'));
tb = poly2sym([cvt(1),cvt(2),cvt(3),cvt(4),cvt(5)],t);
nb = poly2sym([cvn(1),cvn(2),cvn(3),cvn(4),cvn(5),cvn(6),cvn(7),cvn(8)],t);
bb = poly2sym([cvb(1),cvb(2),cvb(3)],t);

%---Part 1---
%(a) GS to Mothership Vect
R_GStoMS = [t.*cos(t) t.*sin(2*t) t];

%(b) Mothership to MAV Vect
%Switch these comments for --P7--
R_MStoMAV = [cos(t) sin(2*t) cos(2*t)];
%R_MStoMAV = [tb nb bb];

% Derivatives
dR = diff(R_GStoMS)/diff(t);
ddR = diff(R_GStoMS,2)/diff(t);

dr_bar = diff(R_MStoMAV)/diff(t);
ddr_bar = diff(R_MStoMAV,2)/diff(t);


%---Part 2---
% Frenet Vectors
t1 = dR/norm(dR);
b1 = cross(dR,ddR)/norm(cross(dR,ddR));
n1 = cross(b1,t1);

t2 = dr_bar/norm(dr_bar);
b2 = cross(dr_bar,ddr_bar)/norm(cross(dr_bar,ddr_bar));
n2 = cross(b2,t2);

% Prep Plots
x1 = subs(R_GStoMS(1),t,time);
y1 = subs(R_GStoMS(2),t,time);
z1 = subs(R_GStoMS(3),t,time);

x2 = subs(R_MStoMAV(1),t,time);
y2 = subs(R_MStoMAV(2),t,time);
z2 = subs(R_MStoMAV(3),t,time);

t_GStoMS = subs(t1,t,time);
n_GStoMS = subs(n1,t,time);
b_GStoMS = subs(b1,t,time);

t_MStoMAV = subs(t2,t,time);
n_MStoMAV = subs(n2,t,time);
b_MStoMAV = subs(b2,t,time);

% Requested Plots
figure
plot3(x1,y1,z1,'color','y'),hold on
title('Groundstation to Mothership')
quiver3(x1,y1,z1,n_GStoMS(:,1),n_GStoMS(:,2),n_GStoMS(:,3),'color','b'),hold on
quiver3(x1,y1,z1,b_GStoMS(:,1),b_GStoMS(:,2),b_GStoMS(:,3),'color','r'),hold on
quiver3(x1,y1,z1,t_GStoMS(:,1),t_GStoMS(:,2),t_GStoMS(:,3),'color','g'),hold on
grid,daspect([1 1 1]),axis vis3d
legend('Curve','Normal','Binormal','Tangent')

figure
plot3(x2,y2,z2,'color','y'),hold on
title('Mothership to MAV')
quiver3(x2,y2,z2,n_MStoMAV(:,1),n_MStoMAV(:,2),n_MStoMAV(:,3),'color','b'),hold on
quiver3(x2,y2,z2,b_MStoMAV(:,1),b_MStoMAV(:,2),b_MStoMAV(:,3),'color','r'),hold on
quiver3(x2,y2,z2,t_MStoMAV(:,1),t_MStoMAV(:,2),t_MStoMAV(:,3),'color','g'),hold on
grid,daspect([1 1 1]),axis vis3d
legend('Curve','Normal','Binormal','Tangent')


%---Part 3---
speed1 = norm(dR);
speed2 = norm(dr_bar);

at1 = dot(ddR,dR)/norm(dR);
at2 = dot(ddr_bar,dr_bar)/norm(dr_bar);

an1 = norm(cross(dR,ddR))/norm(dR);
an2 = norm(cross(dr_bar,ddr_bar))/norm(dr_bar);

s_GStoMS = subs(speed1,t,time);
s_MStoMAV = subs(speed2,t,time);
at_GStoMS = subs(at1,t,time);
at_MStoMAV = subs(at2,t,time);
an_GStoMS = subs(an1,t,time);
an_MStoMAV = subs(an2,t,time);

figure
subplot(2,3,1)
plot(time,s_GStoMS)
xlim([0,2*pi])
title('Speed GS to MS')
subplot(2,3,2)
plot(time,at_GStoMS)
xlim([0,2*pi])
title('Tangential Acc GS to MS')
subplot(2,3,3)
plot(time,an_GStoMS)
xlim([0,2*pi])
title('Normal Acc GS to MS')

subplot(2,3,4)
plot(time,s_MStoMAV,'color','r')
xlim([0,2*pi])
title('Speed MS to MAV')
subplot(2,3,5)
plot(time,at_MStoMAV,'color','r')
xlim([0,2*pi])
title('Tangential Acc MS to MAV')
subplot(2,3,6)
plot(time,an_MStoMAV,'color','r')
xlim([0,2*pi])
title('Normal Acc MS to MAV')


%---Part 4---
% Transformation Vector
G = [t1(1) t1(2) t1(3);
     n1(1) n1(2) n1(3);
     b1(1) b1(2) b1(3)];

% GS to MAV Vector
r = R_GStoMS + (G*R_MStoMAV')';

% Euler Angles
phi = atan(G(1,3)/-G(2,3));
theta = atan(sqrt(1-(G(3,3)^2))/G(3,3));
psi = atan(G(3,1)/G(3,2));

phid = diff(phi)/diff(t);
thetad = diff(theta)/diff(t);
psid = diff(psi)/diff(t);

wbf = [phid*sin(theta)*sin(psi)+thetad*cos(psi) phid*sin(theta)*cos(psi)-thetad*sin(psi) (phid*cos(theta)+psid)];
win = [psid*sin(theta)*sin(phi)+thetad*cos(phi) -psid*sin(theta)*cos(phi)+thetad*sin(phi) (psid*cos(theta)+phid)];

w = wbf;
alpha = diff(w)/diff(t);


%---Part 5---
v = dR + (G*dr_bar')' + cross(w,(G*R_MStoMAV')');

a = ddR + (G*ddr_bar')' + cross(alpha,r) + cross(2*w,(G*dr_bar')') + cross(w,cross(w,r));


%---Part 6---
dr = diff(r)/diff(t);

ddr = diff(dr)/diff(t);

% Part 5&6 Comparison Plots
rVect = feval(matlabFunction(r),time);
rVelo = feval(matlabFunction(dr),time);
rdot = feval(matlabFunction(v),time);

% Speed, Normal&Tangential Acceleration Comparisons
speed3 = norm(dr);
speed4 = norm(v);
rAcc = feval(matlabFunction(ddr),time);
rdubdot = feval(matlabFunction(a),time);
at3 = dot(ddr,dr)/norm(dr);
at4 = dot(a,v)/norm(v);
an3 = norm(cross(dr,ddr))/norm(dr);
an4 = norm(cross(v,a))/norm(v);

% Plot Prep
s_GStoMAV1 = feval(matlabFunction(speed3),time);
s_GStoMAV2 = feval(matlabFunction(speed4),time);
at_GStoMAV1 =feval(matlabFunction(at3),time);
at_GStoMAV2 = feval(matlabFunction(at4),time);
an_GStoMAV1 = feval(matlabFunction(an3),time);
an_GStoMAV2 = feval(matlabFunction(an4),time);

% Plot
hold on
figure
subplot(2,3,1)
plot(time,s_GStoMAV1)
xlim([0,2*pi])
title('P6. Speed GS to MAV')
subplot(2,3,2)
plot(time,at_GStoMAV1)
xlim([0,2*pi])
title('P6. Tangential Acc GS to MAV')
subplot(2,3,3)
plot(time,an_GStoMAV1)
xlim([0,2*pi])
title('P6. Normal Acc GS to MAV')

subplot(2,3,4)
plot(time,s_GStoMAV2,'color','r')
xlim([0,2*pi])
title('P5. Speed GS to MAV')
subplot(2,3,5)
plot(time,at_GStoMAV2,'color','r')
xlim([0,2*pi])
title('P5. Tangential Acc GS to MAV')
subplot(2,3,6)
plot(time,an_GStoMAV2,'color','r')
xlim([0,2*pi])
title('P5. Normal Acc GS to MAV')

figure
hold on
plot3(rVect(:,1), rVect(:,2), rVect(:,3))
hold on 
plot3(x1, y1, z1)
title('Groundstation to MAV')

figure
hold on
plot3(rVelo(:,1), rVelo(:,2), rVelo(:,3))
hold on 
plot3(rdot(:,1), rdot(:,2), rdot(:,3))
legend('part 5 velocity','part 6 velocity')
title('MAV Velocity from GS')

figure
hold on
plot3(rAcc(:,1), rAcc(:,2), rAcc(:,3))
hold on 
plot3(rdubdot(:,1), rdubdot(:,2), rdubdot(:,3))
legend('part 5 acc','part 6 acc')
title('MAV Acceleration from GS')

%{
rVect = double(subs(r,t,time));
rVelo = double(subs(dr,t,time));
rAcc = double(subs(ddr,t,time));

rdot = double(subs(v,t,time));
rdubdot = double(subs(a,t,time));

figure
hold on
plot3(rVect(:,1), rVect(:,2), rVect(:,3))
hold on 
plot3(x1, y1, z1)
title('Groundstation to MAV')

figure
hold on
plot3(rVelo(:,1), rVelo(:,2), rVelo(:,3))
hold on 
plot3(rdot(:,1), rdot(:,2), rdot(:,3))
legend('part 5 velocity','part 6 velocity')
title('MAV Velocity from GS')

figure
hold on
plot3(rAcc(:,1), rAcc(:,2), rAcc(:,3))
hold on 
plot3(rdubdot(:,1), rdubdot(:,2), rdubdot(:,3))
legend('part 5 acc','part 6 acc')
title('MAV Acceleration from GS')
%}
toc