%
% test hphc diff (~1e-14)
% 2022/11/13
% seems theta_s,phi_s,theta_k,phi_k value are a little different 



%  read LDC hdf5 file
Author = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/Author');
Model = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/Approximant');

% Cadence
dt = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/Cadence');
%Cadence = hdf5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/Cadence');
% ObservationDuration
dur= h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/ObservationDuration');
% PolarAngleOfSpin
theta_k = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/PolarAngleOfSpin');
%AzimuthalAngleOfSpin
phi_k = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/AzimuthalAngleOfSpin');
%Distance
D = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/Distance');
%EclipticLatitude
theta_s = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/EclipticLatitude');
%EclipticLongitude
phi_s = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/EclipticLongitude');
%InitialAlphaAngle
alp0 = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/InitialAlphaAngle');
%InitialAzimuthalOrbitalFrequency
nu0 = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/InitialAzimuthalOrbitalFrequency');
%InitialAzimuthalOrbitalPhase
phi0 = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/InitialAzimuthalOrbitalPhase');
%InitialEccentricity
e0 = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/InitialEccentricity');
%InitialTildeGamma
gam0 = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/InitialTildeGamma');
%LambdaAngle
lamd = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/LambdaAngle');
%MassOfCompactObject
u = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/MassOfCompactObject');
%MassOfSMBH
M = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/MassOfSMBH');
%PlungeTime
t_pl= h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/PlungeTime');
%SMBHspin
aa = h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/SMBHspin');
% t, h+, hx
%hphcData
hphcData =  h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/hphcData');
% 't, value'
%TDIdata
TDIdata  =  h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/PreProcess/TDIdata');

LDCparams = [u,M,lamd,aa, e0, nu0,theta_s,phi_s,theta_k,phi_k,  phi0,gam0,alp0,  D];
for i=1:14
    fprintf("%d    %.18e\n", i, LDCparams(i));
end

% 2022/11/13 read LDC params as %.18e
% 1    2.948999954776579813e+01
% 2    1.134944869275097968e+06
% 3    2.142199999999910620e+00
% 4    9.697000000000000064e-01
% 5    4.989444892403920306e-01
% 6    2.232796975919999927e+00
% 7    1.522100180545765014e+00
% 8    3.946698166040000011e+00
% 9    2.041531919926725891e+00
% 10    7.380463140752471611e-04
% 11    5.659686197804838059e+00
% 12    2.286566522026621529e-01
% 13    1.175479035975818931e+00
% 14    5.235888314207546301e+00