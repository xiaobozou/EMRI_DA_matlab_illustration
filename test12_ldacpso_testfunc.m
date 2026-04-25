%% Test script for PSO for least square

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
Dt = D*Gpct;

dt_ode = 15360;
dt_orbit = 61440; %　61440 = 2^12*15

tVecl = (0:N_fft)*dt;
Nl = length(tVecl);
tVecs_lisaorbit = 0:dt_orbit:dur;
tVecs_ode = 0:dt_ode:dur;

% -----------------  ODE->HAK->TDI -----------------%
do_shift = 1; % shift template XYZ, start from 42

%% injected params
u                       =29.489999547765798;
M                       = 1.134944869275098e+06;
lam                     =2.142199999999911;
aa                      =0.969700000000000;
e0                      =0.228656652202662;
nu0                     =7.380463140752472e-04;

theta_s                 =0.498944489240392;
phi_s                   =2.232796975920000;
theta_k                 =1.522100180545765;
phi_k                   =3.946698166040000;

phi0                    =2.041531919926726;
gam0                    =5.659686197804838;
alp0                    =1.175479035975819;

p_true = [u,M,lam,aa,e0,nu0,theta_s,phi_s,theta_k,phi_k,phi0,gam0,alp0];
%% Run PSO
% The fitness function called is LINCHIRPMFFITFUNC. The FFT of the data
% realization is passed through the fitness function input parameter
% structure (2nd input argument). Make sure that the search range includes
% the signal parameters as specified in the generation of the data
% realization. The first search parameter is f0 and the second is f1.
nDim = 13;
f1 = 1-1e-3;
f2 = 1+1e-3;
ffparams = struct('rmin',[f1*p_true],...
                  'rmax',[f2*p_true],...
                  'p_true',p_true,...
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
%%
% Fitness function handle.
fitFuncHandle = @(x) ldacpsotestfunc(x,ffparams);
%%
% Call PSO.
psoOut = ldacpso(fitFuncHandle,nDim);

%% Estimated parameters
% Best standardized and real coordinates found.
stdCoord = psoOut.bestLocation;
bestSNR  = psoOut.bestSNR;
PSO_dump = psoOut.PSO_dump;

fprintf('1 bestfit y=%f\n',bestSNR);
[fitval,realCoord] = fitFuncHandle(stdCoord);
fprintf('2 bestfit y=%f\n',fitval);
for id = 1:13
    fprintf('%d  Estimated=%f; Real =%f\n',id,realCoord(id),p_true(id));
end


% f1 = 1-1e-10;
% f2 = 1+1e-10;
% 4 particle +4 step
% ffparams = struct('rmin',[f1*p_true],...
%                   'rmax',[f2*p_true],...
% 2 bestfit y=-47.866671
% 1  Estimated=29.490000; Real =29.490000
% 2  Estimated=1134944.869283; Real =1134944.869275
% 3  Estimated=2.142200; Real =2.142200
% 4  Estimated=0.969700; Real =0.969700
% 5  Estimated=0.228657; Real =0.228657
% 6  Estimated=0.000738; Real =0.000738
% 7  Estimated=0.498944; Real =0.498944
% 8  Estimated=2.232797; Real =2.232797
% 9  Estimated=1.522100; Real =1.522100
% 10  Estimated=3.946698; Real =3.946698
% 11  Estimated=2.041532; Real =2.041532
% 12  Estimated=5.659686; Real =5.659686
% 13  Estimated=1.175479; Real =1.175479


% f1 = 1-1e-3;
% f2 = 1+1e-3;
% 40 particle +100 step
% 2 bestfit y=-15.557525
% 1  Estimated=29.497364; Real =29.490000
% 2  Estimated=1134144.190600; Real =1134944.869275
% 3  Estimated=2.143815; Real =2.142200
% 4  Estimated=0.970365; Real =0.969700
% 5  Estimated=0.228759; Real =0.228657
% 6  Estimated=0.000738; Real =0.000738
% 7  Estimated=0.498590; Real =0.498944
% 8  Estimated=2.233716; Real =2.232797
% 9  Estimated=1.521836; Real =1.522100
% 10  Estimated=3.945962; Real =3.946698
% 11  Estimated=2.042600; Real =2.041532
% 12  Estimated=5.664230; Real =5.659686
% 13  Estimated=1.174402; Real =1.175479