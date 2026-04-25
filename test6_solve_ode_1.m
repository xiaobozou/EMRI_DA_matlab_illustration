%
% test ode 
%
clc
clear

num_harmonics = 25;
%% LDC params
u       =     2.948999954776579813e+01;
M       =     1.134944869275097968e+06;
lam     =     2.142199999999910620e+00;
aa      =     9.697000000000000064e-01;
theta_s =     4.989444892403920306e-01;
phi_s   =     2.232796975919999927e+00;
theta_k =     1.522100180545765014e+00;
phi_k   =     3.946698166040000011e+00;
phi0    =     2.041531919926725891e+00;
nu0     =     7.380463140752471611e-04;
gam0    =     5.659686197804838059e+00;
e0      =     2.286566522026621529e-01;
alp0    =     1.175479035975818931e+00;
D       =     5.235888314207546301e+00;

dur                     =2392003*15;
dt = 15;
N_fft = dur/dt;
dt_ode = 15360;
dt_orbit = 61440; %　61440 = 2^12*15

tVecl = (0:N_fft)*dt;
Nl = length(tVecl);
tVecs_lisaorbit = 0:dt_orbit:dur;
tVecs_ode = 0:dt_ode:dur;

% unit transfer
% constant in SI unit
G = 6.67408e-11;   % m^3 kg^-1 s^-2
c = 299792458;       % m
% const double LDC_pc = 3.08567758149136727e+16;   /*from LDCConstant.py*/
% const double LISAWP_PC_SI =       3.0856775807e16; /* Parsec, m *, Constant.cc/
pc = 3.0856775807e16;
kpc = 1e3*pc;
Mpc = 1e6*pc;
Gpc = 1e9*pc;
pct = pc/c;
Mpct = Mpc/c;
Gpct = Gpc/c;
%solarmassintime = solarmass*G/(c*c*c);   %time~sec, mass~sec
% 2022/11/14
% // const double LISAWP_MTSUN_SI =    4.92549095e-6;   /* Geometrized solar mass, s, Constant.cc */
solarmassintime = 4.92549095e-6; 

ut = u*solarmassintime;
Mt = M*solarmassintime;
Dt = D*Gpct;

params1 = struct('eclipticLat',theta_s, ...
         'eclipticLong',phi_s, ...
         'theta_k',theta_k, ...
         'phi_k', phi_k, ...
         'ut',ut, ...
         'Mt',Mt, ...
         'aa', aa, ...
         'lam', lam, ...
         'nu0', nu0, ...
         'phi0', phi0, ...
         'alp0', alp0, ...
         'gam0', gam0 , ...
         'e0', e0, ...
         'Dt', Dt, ...
         'dur', dur, ...
         'dt', dt,...
         'dt_ode',dt_ode,...
         'tVecl', tVecl,...
         'tVecs_ode', tVecs_ode,...
         'tVecs_lisaorbit', tVecs_lisaorbit);

%% load LDC data
% t, h+, hx
%hphcData
hphcData =  h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/hphcData');
% 't, value'
%TDIdata
TDIdata  =  h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/PreProcess/TDIdata');

t   = hphcData(1,1:Nl);
hp  = hphcData(2,1:Nl);
hc  = hphcData(3,1:Nl);
X   = TDIdata(2,1:Nl);
Y   = TDIdata(3,1:Nl);
Z   = TDIdata(4,1:Nl);
%%  fidutial orbit and response 
[xts,yts,zts,n1s,n2s,n3s,R0s,R1s,R2s,R3s,q1s,q2s,q3s,Lt]  = func_LISAOrbit(tVecs_lisaorbit);
[Fps,Fcs,k] = func_FpFc(theta_s,phi_s,tVecs_lisaorbit,n1s,n2s,n3s);
fprintf("%f    %f    %f\n",k(1),k(2),k(3));

%% interp response
[Fp,Fc] =interp_FpFc(tVecs_lisaorbit,Fps,Fcs,tVecl);
[kn,kR] = interp_knkR(k,tVecs_lisaorbit,R1s,R2s,R3s,n1s,n2s,n3s,tVecl);

%% TDI 1.0
slr_X = [1,-3, 2; 2, 3, 1;1,2, 3;3,-2, 1;1,2, 3;3,-2, 1;1,-3, 2;2, 3, 1];
slr_Y = [2,-1, 3; 3, 1, 2;2,3, 1;1,-3, 2;2,3, 1;1,-3, 2;2,-1, 3;3, 1, 2];
slr_Z = [3,-2, 1; 1, 2, 3;3,1, 2;2,-1, 3;3,1, 2;2,-1, 3;3,-2, 1;1, 2, 3];
slr_idx = [slr_X;slr_Y;slr_Z];
shiftL = [3*Lt,2*Lt,1*Lt,0*Lt,3*Lt,2*Lt,1*Lt,0*Lt];
shiftL1 = [shiftL,shiftL,shiftL];

%% ODEsolver: ode45   自适应步长算法
AK_ode = @(t,y)Orbit_ODE(t,y,ut,Mt,lam,aa);
y0 = [phi0, nu0, gam0, e0, alp0];
tspan = [0, dur];
% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9, ...
%                  'Stats', 'on', 'OutputFcn', @odeplot);

% --------------- 1 ------------------%
% 相对容差（RelTol）的默认值是 1e-3（即0.001）
% 绝对容差（AbsTol）的默认值是 1e-6
% sol=ode45(AK_ode, tspan, y0);
% size(sol.y) = 5    18

% --------------- 2 ------------------%
% options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);
% sol=ode45(AK_ode, tspan, y0, options);
% size(sol.y) = 5    22
% --------------- 3 ------------------%
options = odeset('RelTol', 1e-9, 'AbsTol', 1e-12);
sol=ode45(AK_ode, tspan, y0, options);
% size(sol.y) = 5    63


% 插值计算
y_eval = deval(sol, tVecl);

[hp1,hc1]=func_get_HAK(theta_s,phi_s,theta_k,phi_k,ut,Mt,aa,lam,Dt,y_eval(1,:),y_eval(2,:),y_eval(3,:),y_eval(4,:),y_eval(5,:),num_harmonics);
[X1,Y1,Z1] = func_XYZ_Mohanty(hp1,hc1,Fp,Fc,kn,kR,tVecl,Lt,slr_idx,shiftL1);
X1 = circshift(X1, -41);
Y1 = circshift(Y1, -41);
Z1 = circshift(Z1, -41);

figure(1)
subplot(5,2,1)
plot(sol.x, sol.y(1,:),'o')
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('$\phi$','Interpreter', 'latex')
subplot(5,2,2)
plot(tVecl, y_eval(1,:))
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('$\phi$','Interpreter', 'latex')
subplot(5,2,3)
plot(sol.x, sol.y(2,:),'o')
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('$\nu$','Interpreter', 'latex')
subplot(5,2,4)
plot(tVecl, y_eval(2,:))
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('$\nu$','Interpreter', 'latex')
subplot(5,2,5)
plot(sol.x, sol.y(3,:),'o')
xlabel('time(sec)')
ylabel('$\gamma$','Interpreter', 'latex')
set(gca,'FontSize',20)  
subplot(5,2,6)
plot(tVecl, y_eval(3,:))
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('$\gamma$','Interpreter', 'latex')
subplot(5,2,7)
plot(sol.x, sol.y(4,:),'o')
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('$e$','Interpreter', 'latex')
subplot(5,2,8)
plot(tVecl, y_eval(4,:))
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('$e$','Interpreter', 'latex')
subplot(5,2,9)
plot(sol.x, sol.y(5,:),'o')
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('$\alpha$','Interpreter', 'latex')
subplot(5,2,10)
plot(tVecl, y_eval(5,:))
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('$\alpha$','Interpreter', 'latex')
screenSize = get(0, 'ScreenSize'); % Get screen resolution [1,4](@ref)  
set(gcf, 'OuterPosition', screenSize); % Cover entire screen  
sgt = sgtitle('Test ode solver, 2026/04/23','Color','k','FontSize',30,'Interpreter', 'latex');
saveas(gcf,'./Figures/test_6_1.png')

%% ODEsolver: RK  fixed步长算法
[phit,nut,gamt,et,alpt,len] = solve_ode(phi0,nu0,gam0,e0,alp0,tVecl,ut,Mt,lam,aa);
[hp2,hc2]=func_get_HAK(theta_s,phi_s,theta_k,phi_k,ut,Mt,aa,lam,Dt,phit,nut,gamt,et,alpt,num_harmonics);
[X2,Y2,Z2] = func_XYZ_Mohanty(hp2,hc2,Fp,Fc,kn,kR,tVecl,Lt,slr_idx,shiftL1);
X2 = circshift(X2, -41);
Y2 = circshift(Y2, -41);
Z2 = circshift(Z2, -41);
%% corr hphc/XYZ test
fprintf("--------------------------\n");
M = corrcoef(hp,hp1);
fprintf('ODE45, corr hp = %f\n',M(1,2));
M = corrcoef(hc,hc1);
fprintf('ODE45, corr hc = %f\n',M(1,2));
M = corrcoef(X,X1);
fprintf('ODE45, corr X = %f\n',M(1,2));
M = corrcoef(Y,Y1);
fprintf('ODE45, corr Y = %f\n',M(1,2));
M = corrcoef(Z,Z1);
fprintf('ODE45, corr Z = %f\n',M(1,2));
fprintf("--------------------------\n");

M = corrcoef(hp,hp2);
fprintf('RK, corr hp = %f\n',M(1,2));
M = corrcoef(hc,hc2);
fprintf('RK, corr hc = %f\n',M(1,2));
M = corrcoef(X,X2);
fprintf('RK, corr X = %f\n',M(1,2));
M = corrcoef(Y,Y2);
fprintf('RK, corr Y = %f\n',M(1,2));
M = corrcoef(Z,Z2);
fprintf('RK, corr Z = %f\n',M(1,2));
fprintf("--------------------------\n");

figure(2)
subplot(5,3,1)
plot(t,hp)
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('hp','Interpreter', 'latex')
legend('LDC')
subplot(5,3,2)
plot(t,hp-hp1)
set(gca,'FontSize',20) 
xlabel('time(sec)')
ylabel('hp','Interpreter', 'latex')
legend('diff LDC/ODE45')
subplot(5,3,3)
plot(t,hp-hp2)
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('hp','Interpreter', 'latex')
legend('diff LDC/RK')
subplot(5,3,4)
plot(t,hc)
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('hc','Interpreter', 'latex')
legend('LDC')
subplot(5,3,5)
plot(t,hc-hc1)
set(gca,'FontSize',20) 
xlabel('time(sec)')
ylabel('hc','Interpreter', 'latex')
legend('diff LDC/ODE45')
subplot(5,3,6)
plot(t,hc-hc2)
set(gca,'FontSize',20)
xlabel('time(sec)')
ylabel('hc','Interpreter', 'latex')
legend('diff LDC/RK')
subplot(5,3,7)
plot(t,X)
set(gca,'FontSize',20) 
xlabel('time(sec)')
ylabel('X','Interpreter', 'latex')
legend('LDC')
subplot(5,3,8)
plot(t,X-X1)
set(gca,'FontSize',20) 
xlabel('time(sec)')
ylabel('X','Interpreter', 'latex')
legend('diff LDC/ODE45')
subplot(5,3,9)
plot(t,X-X2)
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('X','Interpreter', 'latex')
legend('diff LDC/RK')
subplot(5,3,10)
plot(t,Y)
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('Y','Interpreter', 'latex')
legend('LDC')
subplot(5,3,11)
plot(t,Y-Y1)
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('Y','Interpreter', 'latex')
legend('diff LDC/ODE45')
subplot(5,3,12)
plot(t,Y-Y2)
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('Y','Interpreter', 'latex')
legend('diff LDC/RK')
subplot(5,3,13)
plot(t,Z)
set(gca,'FontSize',20) 
xlabel('time(sec)')
ylabel('Z','Interpreter', 'latex')
legend('LDC')
subplot(5,3,14)
plot(t,Z-Z1)
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('Z','Interpreter', 'latex')
legend('diff LDC/ODE45')
subplot(5,3,15)
plot(t,Z-Z2)
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('Z','Interpreter', 'latex')
legend('diff LDC/RK')
screenSize = get(0, 'ScreenSize'); % Get screen resolution [1,4](@ref)  
set(gcf, 'OuterPosition', screenSize); % Cover entire screen  
sgt = sgtitle('Test ode solver, 2026/04/23','Color','k','FontSize',30,'Interpreter', 'latex');
saveas(gcf,'./Figures/test_6_2.png')
