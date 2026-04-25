% test HAK 14D
% have solved ode already, used as input of func_get_myHAK_4.m
% 2022/01/20

% reduce the loop in harmonics summation
% not store amplitude\&phase of each harmonics
% 2022/02/25 update

% preparing the length N=2^22+1 for hphz,XYZ
% 2022/03/22

clc
clear

%% load params, LDC-1.2
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
% #########################################################
dur                     =62914560;
dt                      =15.0;


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

% 1 yr
%T = 3600*24*365*1.0;  % 1 yr
%T = 2^21*15; % 1 yr

% 2 yr
%T = 3600*24*365*2.0;  % 2 yr
%N = 2^22;
%T = 2^22*15; % 2 yr
dt = 15;
dt_ode = 15360;

N_fft = dur/dt;
tVecl = (0:N_fft)*dt;  % 2^22 +1
% tVecs_ode(end)-T = 0  2022/03/22 check
tVecs_ode = 0:dt_ode:dur;
Ns_ode = length(tVecs_ode);
Nl = length(tVecl);
%% solving ode
%  check_ode_start
start_ok = check_ode_start(e0,nu0,Mt);
fprintf("start_ok = %d\n",start_ok);

% solve ode
if  start_ok == 1   
    % solve ode  in deltat
    phit1  = zeros(1,Ns_ode);
    nut1   = zeros(1,Ns_ode);
    gamt1  = zeros(1,Ns_ode);
    et1    = zeros(1,Ns_ode);
    alpt1  = zeros(1,Ns_ode);
    [phit1,nut1,gamt1,et1,alpt1,len] = solve_ode(phi0,nu0,gam0,e0,alp0,tVecs_ode,ut,Mt,lam,aa);
    fprintf("for deltat, length = %d,  t_pl = %f\n", len, tVecs_ode(len));
    % interp ode to dt
    %-----------------------------------------------------------------
    %--pay attention to the length and number of intervals, 2022/11/04
    %sz_t = ceil((len-1)*deltat/dt);
    %-----------------------------------------------------------------
    len2 = ceil((len-1)*dt_ode/dt) + 1;
    fprintf("for dt, length = %d,  t_pl = %f\n", len2, tVecl(len2)); 
    
    phit = zeros(1,N_fft);
    nut  = zeros(1,N_fft);
    gamt = zeros(1,N_fft);
    et   = zeros(1,N_fft);
    alpt = zeros(1,N_fft);
    [phit(1:len2),nut(1:len2),gamt(1:len2),et(1:len2),alpt(1:len2)] = interp_ode(tVecs_ode(1:len),phit1(1:len),nut1(1:len),gamt1(1:len),et1(1:len),alpt1(1:len),tVecl(1:len2));

    % take interpolated ode to 14D hphc generation
    % 25 harmonics
    num_harmonics = 25;
    % 10 harmonics
    %num_harmonics = 10;
    % 5 harmonics
    %num_harmonics = 5; 
    % 3 harmonics
    %num_harmonics = 3;  
    [hp,hc]=func_get_HAK(theta_s,phi_s,theta_k,phi_k,ut,Mt,aa,lam,Dt,phit,nut,gamt,et,alpt,num_harmonics);
else
    fprintf("can not pass the start_ode_check from the given initial.\n");
end
%% corr hphc test
% t, h+, hx
%hphcData
hphcData =  h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/hphcData');
% 't, value'
%TDIdata
% TDIdata  =  h5read('/home/xiaobo/LDC2_Radler_v2/LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/PreProcess/TDIdata');

t2   = hphcData(1,:);
hp2  = hphcData(2,:);
hc2  = hphcData(3,:);

M = corrcoef(hp2,hp);
fprintf('corr hp = %f\n',M(1,2));

M = corrcoef(hc2,hc);
fprintf('corr hc = %f\n',M(1,2));

%% plot
segment = 100000:105000;
figure
subplot(2,2,1)
plot(tVecl(1:N_fft),hp)
set(gca,'fontsize',20);
xlabel('time(sec)')
ylabel('hp')
legend('full')
subplot(2,2,2)
plot(tVecl(1:N_fft),hc)
set(gca,'fontsize',20);
xlabel('time(sec)')
ylabel('hc')
legend('full')
subplot(2,2,3)
plot(tVecl(segment),hp(segment))
set(gca,'fontsize',20);
xlabel('time(sec)')
ylabel('hp')
legend('segment')
subplot(2,2,4)
plot(tVecl(segment),hc(segment))
set(gca,'fontsize',20);
xlabel('time(sec)')
ylabel('hc')
legend('segment')
screenSize = get(0, 'ScreenSize'); % Get screen resolution [1,4](@ref)  
set(gcf, 'OuterPosition', screenSize); % Cover entire screen  
sgt = sgtitle('Test AK waveform, 2026/04/23','Color','k','FontSize',30,'Interpreter', 'latex');
saveas(gcf,'./Figures/test_1_1.png')