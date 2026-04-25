function [Sx,Sae,St] = LDC_psd(fvec)

c = 299792458.0 ;
%## Global
Soms_d = struct('Proposal',(10.e-12)^2, 'SciRDv1', (15.e-12)^2, 'MRDv1', (10.e-12)^2);  %# m^2/Hz
%### Acceleration
Sa_a = struct('Proposal',(3.e-15)**2, 'SciRDv1', (3.e-15)^2, 'MRDv1', (2.4e-15)^2);  %# m^2/sec^4/Hz









end