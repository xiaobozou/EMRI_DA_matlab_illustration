# 2026/04/25

The LDC-1.2 signal and data can be downloaded from https://lisa-ldc.in2p3.fr/challenge1a 

% t, h+, hx
%hphcData
hphcData =  h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/GWSources/EMRI-0/hphcData');
% 't, value'
%TDIdata
% TDIdata  =  h5read('./LDC1-2_EMRI_v2_noiseless.hdf5','/H5LISA/PreProcess/TDIdata');
% TDIdata  =  h5read('./LDC1-2_EMRI_v2.hdf5','/H5LISA/PreProcess/TDIdata');
