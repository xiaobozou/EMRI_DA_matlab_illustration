%
% 2022/07/25 test HAK SNR and loglike13 in both TD and FD
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


%% orbital ode
%  check_ode_start
start_ok = check_ode_start(e0,nu0,Mt);
fprintf("start_ok = %d\n",start_ok);

if start_ok == 0
    fprintf("check_ode_start   fail \n");
end

%% TDI 1.0
slr_X = [1,-3, 2; 2, 3, 1;1,2, 3;3,-2, 1;1,2, 3;3,-2, 1;1,-3, 2;2, 3, 1];
slr_Y = [2,-1, 3; 3, 1, 2;2,3, 1;1,-3, 2;2,3, 1;1,-3, 2;2,-1, 3;3, 1, 2];
slr_Z = [3,-2, 1; 1, 2, 3;3,1, 2;2,-1, 3;3,1, 2;2,-1, 3;3,-2, 1;1, 2, 3];
slr_idx = [slr_X;slr_Y;slr_Z];

delta = pi/10;
% for id_u = -10:10
%     for id_v = 0:10  

        % id_u = -5;
        % id_v = 0;
        % do_shift=1,  SNR FD  =NaN
        % do_shift=1, LLR TD =NaN


        % id_u = 5;
        % id_v = 0;
        % do_shift=1,  SNR FD  =NaN
        % do_shift=1, LLR TD =NaN

        % id_u = -5;
        % id_v = 10;
        % do_shift=1,  SNR FD  =NaN
        % do_shift=1, LLR TD =NaN


        id_u = 5;
        id_v = 10;
        % do_shift=1,  SNR FD  =NaN
        % do_shift=1, LLR TD =NaN
        % innerprod_hh_A_FD  =NaN
        % innerprod_hh_E_FD  =NaN
        % A(4194265) =NaN
        % E(4194265) =NaN
        % X(4194265) =-4.7663e-38
        % Y(4194265) =NaN
        % Z(4194265) =NaN
        % id_anomay = any(isnan(A));  %1
        % id_anomay = any(isinf(A));  %0
        % id_anomay = any(isnan(E));  %1
        % id_anomay = any(isinf(E));  %0
        % id_anomay = find(A==nan);
        % id_anomay = sum(A==nan);  %0
        % id_anomay = sum(isnan(A)); %1
        

        fprintf("---%d---%d\n",id_u,id_v);
 
        phi_s   = id_u*delta;
        theta_s = id_v*delta;

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
        shiftL = [3*Lt,2*Lt,1*Lt,0*Lt,3*Lt,2*Lt,1*Lt,0*Lt];
        shiftL1 = [shiftL,shiftL,shiftL];
        %% interp response 
        [Fp,Fc] = interp_FpFc(tVecs_lisaorbit,Fps,Fcs,tVecl);
        [kn,kR] = interp_knkR(k,tVecs_lisaorbit,R1s,R2s,R3s,n1s,n2s,n3s,tVecl);



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
        num_harmonics = 25;
        [hp,hc]=func_get_HAK(theta_s,phi_s,theta_k,phi_k,ut,Mt,aa,lam,Dt,phit,nut,gamt,et,alpt,num_harmonics);
        tic;
        [X,Y,Z] = func_XYZ_Mohanty(hp,hc,Fp,Fc,kn,kR,tVecl,Lt,slr_idx,shiftL1);
        time = toc;
        fprintf("total XYZ cost = %f sec\n",time);
    
        % if do_shift == 1
        %     X = circshift(X,-41);
        %     Y = circshift(Y,-41);
        %     Z = circshift(Z,-41);
        % end
    
        %% LLR
        [A,E,~] = AET(X(1:N_fft),Y(1:N_fft),Z(1:N_fft));

        % -----------------------------------%
        % 2026/04/25  tracing the LLR NAN anomaly in (thS,phS) plane
        id_anomay = isnan(A);find(id_anomay==1)  %1
        id_anomay = isnan(E);find(id_anomay==1)  %1

        id_anomay = isnan(Z);find(id_anomay==1)  %1
        id_anomay = isnan(Y);find(id_anomay==1)  %1
        id_anomay = isnan(X);find(id_anomay==1)  %0
        fprintf("A(1)=%e\n",A(1));
        fprintf("E(1)=%e\n",E(1));
        fprintf("X(1)=%e\n",X(1));
        fprintf("Y(1)=%e\n",Y(1));
        fprintf("Z(1)=%e\n",Z(1));
        % -----------------------------------%

        Af  = fft(A)*dt;   Ef = fft(E)*dt; 
        innerprod_hh_A_FD = innerproduct_fre_2(Af,Af,Sae2,df,id_cut);
        innerprod_hh_E_FD = innerproduct_fre_2(Ef,Ef,Sae2,df,id_cut);
        SNR2_fd = innerprod_hh_A_FD+innerprod_hh_A_FD;
        fprintf("do_shift=%d,  SNR FD  =%f\n", do_shift, sqrt(SNR2_fd));
        
        innerprod_dh_A_TD = 2*sum(A.*Adw2)*dt;
        innerprod_dh_E_TD = 2*sum(E.*Edw2)*dt;
        LLR2_TD =  (innerprod_dh_A_TD+innerprod_dh_A_TD)^2/(innerprod_hh_A_FD+innerprod_hh_A_FD);
        fprintf("do_shift=%d, LLR TD =%f\n",do_shift, sqrt(LLR2_TD));

        LLR_2D = sqrt(LLR2_TD);