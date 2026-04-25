%
% test snr of pure noise, that is the snr error estiamtion
%  2022/11/17
%
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
model = 'mldc';
%model = 'Proposal';
%model = 'SciRDv1'; 
%model = 'MRDv1';
%fprintf("noise model is %s\n", model);

unit = 'relativeFrequency';
%freqVec = fftfreqvec1(N_fft,fs);
[Sx2,Sae2,St2] = noisepsd_2(N_fft,fs,model,unit);
[Sx,Sae,St] = noisepsd_1(fVec,model,unit);



%% pure  noise
TDIdata  =  h5read('./LDC1-2_EMRI_v2.hdf5','/H5LISA/PreProcess/TDIdata');
%t     = TDIdata(1,1:end-1);
X1_LC2 = TDIdata(2,1:end-1);
Y1_LC2 = TDIdata(3,1:end-1);
Z1_LC2 = TDIdata(4,1:end-1);

TDIdata  =  h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/PreProcess/TDIdata');
X2_LC2 = TDIdata(2,1:end-1);
Y2_LC2 = TDIdata(3,1:end-1);
Z2_LC2 = TDIdata(4,1:end-1);


% pure noise XYZ
X = X1_LC2 - X2_LC2;
Y = Y1_LC2 - Y2_LC2;
Z = Z1_LC2 - Z2_LC2;

%% snr
snr_all = func_snr(dt, X, Y, Z,df,  Sae2, id_cut);
fprintf("pure noise, snr=%f\n",snr_all);

%pure noise, snr=5340.257326
% meaningless





