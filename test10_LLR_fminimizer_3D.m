%
% 2022/07/25 test fminimizer over 3D initial phase
% interp Fp,Fc,kn,kR
% interp ode


clc
clear
%% psd
dur                     =62914560;
dt                      =15.0;
N_fft = dur/dt;
fs = 1.0/dt;
df = fs/N_fft;
fVec = 0:df:fs/2;
fmin = 1e-4;
id_cut = min(find(fVec - fmin > 0));
id_cut

addpath("./LDC_PSD_matlab_1025");
% psd
%model = 'mldc';
%model = 'Proposal';
model = 'SciRDv1'; 
%model = 'MRDv1';
%fprintf("noise model is %s\n", model);

unit = 'relativeFrequency';
%freqVec = fftfreqvec1(N_fft,fs);
[Sx2,Sae2,St2] = noisepsd_2(N_fft,fs,model,unit);
[Sx,Sae,St] = noisepsd_1(fVec,model,unit);


%% load LDC TDI data
%TDIdata  =  h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/PreProcess/TDIdata');
TDIdata  =  h5read('./LDC1-2_EMRI_v2.hdf5','/H5LISA/PreProcess/TDIdata');
%t     = TDIdata(1,1:end-1);
X_LC2 = TDIdata(2,1:N_fft);
Y_LC2 = TDIdata(3,1:N_fft);
Z_LC2 = TDIdata(4,1:N_fft);
clear TDIdata;
[Ad,Ed,Td] = AET(X_LC2,Y_LC2,Z_LC2);
Adf = fft(Ad*dt); Edf = fft(Ed*dt); 
% double whiten the LDC data
temp = fft(Ad);
temp(1) = 0;
temp(2:id_cut) = 0;
temp((end-(id_cut-2)):end) = 0;
Adw2 = ifft(temp./Sae2);

temp = fft(Ed);
temp(1) = 0;
temp(2:id_cut) = 0;
temp((end-(id_cut-2)):end) = 0;
Edw2 = ifft(temp./Sae2);

%% LDC params
theta_s                 =0.498944489240392;
phi_s                   =2.232796975920000;
theta_k                 =1.522100180545765;
phi_k                  =3.946698166040000;
lam                     =2.142199999999911;
aa                      =0.969700000000000;
u                       =29.489999547765798;
M                       = 1.134944869275098e+06;
e0                      =0.228656652202662;
nu0                     =7.380463140752472e-04;
phi0                    =2.041531919926726;
gam0                   =5.659686197804838;
alp0                   =1.175479035975819;
D                       =5.235888314207546;

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

dt_ode = 15360;
dt_orbit = 61440; %　61440 = 2^12*15

tVecl = (0:N_fft)*dt;
Nl = length(tVecl);
tVecs_lisaorbit = 0:dt_orbit:dur;
tVecs_ode = 0:dt_ode:dur;

% -----------------  ODE->HAK->TDI -----------------%
do_shift = 1; % shift template XYZ, start from 42
%do_shift = 0; 

params = struct('theta_s',theta_s, ...
 'phi_s',phi_s, ...
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
 'tVecs_lisaorbit', tVecs_lisaorbit,...
 'Sae2',Sae2,...
 'Adw2',Adw2,...
 'Edw2',Edw2,...
 'id_cut',id_cut);

% function LLR_TD = fitness_10D(x,params)
         % 
         % f = @(x,c) x(1).^2+c.*x(2).^2;  % The parameterized function.
         % c = 1.5;                        % The parameter.
         % X = fminsearch(@(x) f(x,c),[0.3;1])
         % 
options = optimset('PlotFcns', {@optimplotx, ...        % 参数变化
                                    @optimplotfval, ...     % 函数值
                                    @optimplotfunccount}, ... % 函数调用次数
                   'Display', 'iter', ...             
                   'MaxIter', 100, ...         % 最大迭代次数
                   'TolFun', 1e-6, ...         % 函数值容差
                   'TolX', 1e-6, ...           % 变量容差
               'FunValCheck', 'on');       % 检查函数值
fitfuc = @(x) -fitness_10D(x,params);    
x0 = [0,0,0];
[x_opt, fval] = fminsearch(fitfuc, x0, options);
% x_opt =
% 
% %     2.0002   -2.1010    1.1190
% 
% fval =
% 
%   -47.9166