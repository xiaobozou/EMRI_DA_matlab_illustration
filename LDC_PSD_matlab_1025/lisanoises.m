function [Spm, Sop]=lisanoises(f,model,unit)

clight = 299792458.0;    %# Speed of light in vacuum [m.s^-1] (CODATA 2014)
c = clight;



lisaL = 2.5e9;% # LISA's arm meters
lisaLT = lisaL/clight;% # LISA's armn in sec
lisaD = 0.3;%  # TODO check it
lisaP = 2.0;%  # TODO check it

%## Global
Soms_d = struct('Proposal',(10.e-12)^2, 'SciRDv1', (15.e-12)^2, 'MRDv1', (10.e-12)^2);  %# m^2/Hz
%### Acceleration
Sa_a = struct('Proposal',(3.e-15)^2, 'SciRDv1', (3.e-15)^2, 'MRDv1', (2.4e-15)^2);  %# m^2/sec^4/Hz


%if  model == 'mldc'
if strcmp(model,'mldc')
    Spm = 2.5e-48 * (1.0 + (f/1.0e-4).^-2) .* f.^(-2);
    defaultL = 16.6782;
    Sop = 1.8e-37 * (lisaLT/defaultL)^2 * f.^2;
    
    Sa_d = Spm;  %no meaning just for code run.
    Soms_d = Sop;  %no meaning just for code run.
elseif strcmp(model, 'newdrs')  %# lisalight, to be used with lisaL = 1Gm, lisaP = 2
    Spm = 6.00314e-48 * f.^(-2);                                 %# 4.6e-15 m/s^2/sqrt(Hz)
    defaultL = 16.6782;
    defaultD = 0.4;
    defaultP = 1.0;
    Sops = 6.15e-38 * (lisaLT/defaultL)^2 * (defaultD/lisaD)^4 * (defaultP/lisaP);      %# 11.83 pm/sqrt(Hz)
    Sopo = 2.81e-38;                                                                                                         %# 8 pm/sqrt(Hz)
    Sop = (Sops + Sopo) * f.^2;

    Sa_d = Spm;  %no meaning just for code run.
    Soms_d = Sop;  %no meaning just for code run.
elseif strcmp(model,'LCESAcall')
    frq = f
    %### Acceleration noise
    %## In acceleration
    Sa_a1 = Sa_a.Proposal *(1.0 +(0.4e-3./frq).^2+(frq/9.3e-3).^4);
    %## In displacement
    Sa_d = Sa_a1.*(2.*pi*frq).^(-4.);
    %## In relative frequency unit
    Sa_nu = Sa_d.*(2.0*pi*frq/clight).^2;
    Spm =  Sa_nu;

    %### Optical Metrology System
    %## In displacement
    Soms_d1 = Soms_d.Proposal * (1. + (2.e-3./f).^4);
    %## In relative frequency unit
    Soms_nu = Soms_d1.*(2.0*pi*frq/clight).^2;
    Sop =  Soms_nu;
elseif strcmp(model,'Proposal')
    frq = f;
    %### Acceleration noise
    %## In acceleration
    Sa_a1 = Sa_a.Proposal *(1.0 +(0.4e-3./frq).^2).*(1.0+(frq/8e-3).^4);
    %## In displacement
    Sa_d = Sa_a1.*(2.*pi*frq).^(-4.);
    %## In relative frequency unit
    Sa_nu = Sa_d.*(2.0*pi*frq/clight).^2;
    Spm =  Sa_nu;

    %### Optical Metrology System
    %## In displacement
    Soms_d1 = Soms_d.Proposal * (1. + (2.e-3./f).^4);
    %## In relative frequency unit
    Soms_nu = Soms_d1.*(2.0*pi*frq/clight).^2;
    Sop =  Soms_nu;
elseif   strcmp(model,'SciRDv1') || strcmp(model,'MRDv1')
    frq = f;
    %### Acceleration noise
    %## In acceleration
    Sa_a1 = Sa_a.SciRDv1 *(1.0 +(0.4e-3./frq).^2).*(1.0+(frq/8e-3).^4);
    %## In displacement
    Sa_d = Sa_a1.*(2.*pi*frq).^(-4.);
    %## In relative frequency unit
    Sa_nu = Sa_d.*(2.0*pi*frq/clight).^2;
    Spm =  Sa_nu;

    %### Optical Metrology System
    %## In displacement
    Soms_d1 = Soms_d.SciRDv1 * (1. + (2.e-3./frq).^4);
    %## In relative frequency unit
    Soms_nu = Soms_d1.*(2.0*pi*frq/clight).^2;
    Sop =  Soms_nu;
end


if strcmp(unit,'relativeFrequency') 
        Spm = Spm;
        Sop = Sop;
elseif strcmp(unit,'displacement')
       Spm = Sa_d;
       Sop = Soms_d1;
end


end