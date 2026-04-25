% generate XYZ
% interp ode+ interp FpFc, 
% 2022/01/20
clc
clear
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
%T = 2^22*15; % 2 yr


%---------------------------------------
% dt_ode = 15;
%---------------------------------------
% 2022/0105
% bug  for interp ode calling func_get_myHAK_0105.m
% corr A,E, 0.9158  0.9161
%
% when dt_ode = dt, the ode have no interp and should have same result as ode dt.
% why are they different now?
% works well for no interp ode calling func_get_myHAK_3.m
% % corr A,E, 0.9813, 0.9791
%
% when dt_ode = dt,
% works well for for c, see 1225_4 of c
%
% % % # dt_ode = dt
% % % # in matlab
% % % # FIXME bug still exist
% % % # FIXME1  in file func_get_myHAK_0105.m, using wrong function name func_get_myHAK 
% % % # FIXME bug still exist
        %
        %  error  nXp = n*Xp(n)
        %
        %
        % 0105 bug
        % nXp = n*Xp(n+1);
        %
        %
        % shoule be
        %nXp = n*Xp(n+1,:);
        %
        %
% % % #  bug FIXED

%---------------------------------------

%---------------------------------------
%dt_ode = 15360;
%---------------------------------------
%dt_ode = 15;
 % 2022/01/18
 % testing Mohanty's AE snr , c 61.856858, matlab	61.440524
 % result from the significant AE residual for start and end transient.
 % the start transient comes from interp over out of range for time delay.
 % the end transient comes from accumulated inaccuray bring by interp ode.
 % see 1001_Mohanty/0118_ode_calibration
 
 % dt_ode = dt; no ode interp inaccuracy, expexting snr c&matlab come
 % closer.
 % still 61.8 c, 61.4 matlab, tell ode interp error can be ignored.
%---------------------------------------
dt = 15;
dt_ode = 15360;
%dt_ode = 15;
N_fft = dur/dt;
tVecl = (0:N_fft)*dt;
Nl = length(tVecl);

% 2022/03/22
% log2(86400*14) ~ 20.2061 
dt_orbit = 61440;  % 2^12*15 = 61440
tVecs_lisaorbit = 0:dt_orbit:dur;
% tVecs_lisaorbit(end) -T = 0
tVecs_ode = 0:dt_ode:dur;
% tVecs_ode(end)-T = 0

% 13D params
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


%%  fidutial orbit and response 
[xts,yts,zts,n1s,n2s,n3s,R0s,R1s,R2s,R3s,q1s,q2s,q3s,Lt]  = func_LISAOrbit(tVecs_lisaorbit);
[Fps,Fcs,k] = func_FpFc(theta_s,phi_s,tVecs_lisaorbit,n1s,n2s,n3s);
fprintf("%f    %f    %f\n",k(1),k(2),k(3));

%% interp response
[Fp,Fc] =interp_FpFc(tVecs_lisaorbit,Fps,Fcs,tVecl);
[kn,kR] = interp_knkR(k,tVecs_lisaorbit,R1s,R2s,R3s,n1s,n2s,n3s,tVecl);

% -----------------  ODE->HAK->TDI -----------------%
%% orbital ode
%  check_ode_start
start_ok = check_ode_start(e0,nu0,Mt);
fprintf("start_ok = %d\n",start_ok);

if start_ok == 0
    fprintf("check_ode_start   fail \n");

elseif start_ok == 1
        
    % solve ode  in dt_ode
    [phit1,nut1,gamt1,et1,alpt1,len] = solve_ode(phi0,nu0,gam0,e0,alp0,tVecs_ode,ut,Mt,lam,aa);
    fprintf("for dt_ode, len = %d,  tVecs_ode(len) = %f\n", len, tVecs_ode(len));


    % interp ode
    sz_t = (len-1)*dt_ode/dt+1;
    fprintf("for dt, sz_t = %d,  tVecl(sz_t) = %f\n",  sz_t, tVecl(sz_t));
    
    phit = zeros(1,Nl);
    nut  = zeros(1,Nl);
    gamt = zeros(1,Nl);
    et   = zeros(1,Nl);
    alpt = zeros(1,Nl);
    [phit,nut,gamt,et,alpt] = interp_ode(tVecs_ode, phit1,nut1,gamt1,et1,alpt1,tVecl);
   
    %% HAK
    % take interolated ode tp 14D hphc generation
    % 25 harmonics
    num_harmonics = 25;
    % 10 harmonics
    %num_harmonics = 10;
    % 5 harmonics
    %num_harmonics = 5; 
    % 3 harmonics
    %num_harmonics = 3;  
    [hp,hc]=func_get_HAK(theta_s,phi_s,theta_k,phi_k,ut,Mt,aa,lam,Dt,phit,nut,gamt,et,alpt,num_harmonics);

    %% TDI 1.0
    slr_X = [1,-3, 2; 2, 3, 1;1,2, 3;3,-2, 1;1,2, 3;3,-2, 1;1,-3, 2;2, 3, 1];
    slr_Y = [2,-1, 3; 3, 1, 2;2,3, 1;1,-3, 2;2,3, 1;1,-3, 2;2,-1, 3;3, 1, 2];
    slr_Z = [3,-2, 1; 1, 2, 3;3,1, 2;2,-1, 3;3,1, 2;2,-1, 3;3,-2, 1;1, 2, 3];
    slr_idx = [slr_X;slr_Y;slr_Z];
    shiftL = [3*Lt,2*Lt,1*Lt,0*Lt,3*Lt,2*Lt,1*Lt,0*Lt];
    shiftL1 = [shiftL,shiftL,shiftL];

    tic;
    [X,Y,Z] = func_XYZ_Mohanty(hp,hc,Fp,Fc,kn,kR,tVecl,Lt,slr_idx,shiftL1);
    time = toc;
    fprintf("total XYZ cost = %f sec\n",time);

end

%% corr check
% LDC TDI
TDIdata  =  h5read('/home/xiaobozou/Projects_LISA/Data/LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/PreProcess/TDIdata');
t1 = TDIdata(1,:);
X1 = TDIdata(2,:);
Y1 = TDIdata(3,:);
Z1 = TDIdata(4,:);


[A,E,T] = AET(X,Y,Z);
[A1,E1,T1] = AET(X1,Y1,Z1);


% check shift
[c1,lags1] = xcorr(X,X1,1000,'normalized');
[c2,lags2] = xcorr(Y,Y1,1000,'normalized');
[c3,lags3] = xcorr(Z,Z1,1000,'normalized');
[c4,lags4] = xcorr(A,A1,1000,'normalized');
[c5,lags5] = xcorr(E,E1,1000,'normalized');

figure(1)
subplot(1,3,1)
stem(lags1,c1)
set(gca,'FontSize',20)  
subplot(1,3,2)
stem(lags2,c2)
set(gca,'FontSize',20)  
subplot(1,3,3)
stem(lags3,c3)
set(gca,'FontSize',20)  
screenSize = get(0, 'ScreenSize'); % Get screen resolution [1,4](@ref)  
set(gcf, 'OuterPosition', screenSize); % Cover entire screen  
sgt = sgtitle('Test TDI comparing with LDC-1.2, 2026/04/23','Color','k','FontSize',30,'Interpreter', 'latex');
saveas(gcf,'./Figures/test_3_1.png')


[max_c1,id1] = max(c1); 
[max_c2,id2] = max(c2); 
[max_c3,id3] = max(c3);
[max_c4,id4] = max(c4);
[max_c5,id5] = max(c5);
fprintf('-----------------------------------------\n')
fprintf('---X is template ---\n')
fprintf('---X1 is LDC signal ---\n')
fprintf("---xcorr(X,X1) max_corr=%f   lag=%d\n", max_c1,lags1(id1));
fprintf("---xcorr(Y,Y1) max_corr=%f   lag=%d\n", max_c2,lags2(id2));
fprintf("---xcorr(Z,Z1) max_corr=%f   lag=%d\n", max_c3,lags3(id3));
fprintf("---xcorr(A,A1) max_corr=%f   lag=%d\n", max_c4,lags3(id4));
fprintf("---xcorr(E,E1) max_corr=%f   lag=%d\n", max_c5,lags3(id5));
fprintf('-----------------------------------------\n')


%% plot
% --------1, shift XYZ by 615 sec as LC2-----------%
X_sh = circshift(X,-41);
Y_sh = circshift(Y,-41);
Z_sh = circshift(Z,-41);
figure(2)
subplot(3,1,1)
plot(t1,X_sh,t1,X1,t1,X1-X_sh)
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('X')
legend('me','LDC','diff')
subplot(3,1,2)
plot(t1,Y_sh,t1,Y1,t1,Y1-Y_sh)
set(gca,'FontSize',20)  
xlabel('time(sec)')
ylabel('Y')
legend('me','LDC','diff')
subplot(3,1,3)
plot(t1,Z_sh,t1,Z1,t1,Z1-Z_sh)
set(gca,'FontSize',20)
xlabel('time(sec)')
ylabel('Z')
legend('me','LDC','diff')
screenSize = get(0, 'ScreenSize'); % Get screen resolution [1,4](@ref)  
set(gcf, 'OuterPosition', screenSize); % Cover entire screen  
sgt = sgtitle('Test TDI comparing with LDC-1.2, 2026/04/23','Color','k','FontSize',30,'Interpreter', 'latex');
saveas(gcf,'./Figures/test_3_2.png')