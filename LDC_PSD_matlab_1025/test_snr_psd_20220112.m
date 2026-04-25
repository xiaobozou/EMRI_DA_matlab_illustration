% test snr

%clc
%clear

%TDIdata
TDIdata0  =  h5read('/home/xiaobo/LDC2_Radler_v2/LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/PreProcess/TDIdata');
TDIdata1  =  h5read('/home/xiaobo/LDC2_Radler_v2/LDC1-2_EMRI_v2.hdf5','/H5LISA/PreProcess/TDIdata');

X0 = TDIdata0(2,:); Y0 = TDIdata0(3,:); Z0 = TDIdata0(4,:);
[A,E,T] = AET(X0,Y0,Z0);
Nt = length(A);

%id_half = round(Nt/2)+1;
id_half = round(Nt/2);

dt = 15;
fs = 1/dt;
df = fs/Nt;


fmin = 1e-4;  
fmax = 0.1; 


% %% calculating snr using 1-sided psd in frequency domain % (h|h)
% Af = fft_py(A*dt); Ef = fft_py(E*dt); Tf = fft_py(T*dt);
%% calculating snr using 2-sided psd in frequency domain % (h|h)
Af = fft(A*dt); Ef = fft(E*dt); Tf = fft(T*dt);
fvecs_positive = (0:id_half-1)*df;
% half complex fvec
fvec=fftfreqvec1(Nt,fs);
    
id_cut = min(find(fvecs_positive - fmin > 0));
id_cut
% 2-sided PSD
%[Sx,Sae,St] = PSD(fvec);
% for PSD, fvec should be positive
%
% 2022/01/10 test psd
%[Sx,Sae,St] = PSD(abs(fvec));

%model = 'mldc';
model = 'Proposal';
unit = 'relativeFrequency';
%unit = 'displacement';

Sx = noisepsd_X(abs(fvec), model,unit);
Sae = noisepsd_AE(abs(fvec), model,unit);
St = noisepsd_T(abs(fvec), model,unit);

%% all band
y=innerProd_fre_2SidedPSD(Af,Af,Sae,df);
snr_aa = sqrt(y);
y=innerProd_fre_2SidedPSD(Ef,Ef,Sae,df);
snr_ee = sqrt(y);
y=innerProd_fre_2SidedPSD(Tf,Tf,St,df);
snr_tt = sqrt(y);

snr_total =  sqrt(snr_aa^2+snr_ee^2+snr_tt^2);
fprintf("in frequency domain all band\n")
fprintf("using 2-sided PSD\n")
fprintf("orginal signal snr AET =  %f %f %f  %f\n", snr_aa,snr_ee,snr_tt,snr_total);
fprintf("\n");



%% only LISA band

y=innerProd_fre_2SidedPSD(Af(id_cut:end-id_cut),Af(id_cut:end-id_cut),Sae(id_cut:end-id_cut),df);
snr_aa = sqrt(y);
y=innerProd_fre_2SidedPSD(Ef(id_cut:end-id_cut),Ef(id_cut:end-id_cut),Sae(id_cut:end-id_cut),df);
snr_ee = sqrt(y);
y=innerProd_fre_2SidedPSD(Tf(id_cut:end-id_cut),Tf(id_cut:end-id_cut),St(id_cut:end-id_cut),df);
snr_tt = sqrt(y);

snr_total =  sqrt(snr_aa^2+snr_ee^2+snr_tt^2);
fprintf("in frequency domain LISA band\n")
fprintf("using 2-sided PSD\n")
fprintf("orginal signal snr AET =  %f %f %f  %f\n", snr_aa,snr_ee,snr_tt,snr_total);

%
% 2022/01/10 the frequency cut within LISAband change a lot snr_T
%