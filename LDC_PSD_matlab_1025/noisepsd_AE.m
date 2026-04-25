function Sae = noisepsd_AE(fvec,model,unit)


clight = 299792458.0;    %# Speed of light in vacuum [m.s^-1] (CODATA 2014)
c = clight;
lisaL = 2.5e9;% # LISA's arm meters
lisaLT = lisaL/clight;% # LISA's armn in sec

    x = 2.0 * pi * lisaLT * fvec;

    [Spm, Sop] = lisanoises(fvec, model,unit);

    Sae = 8.0 * sin(x).^2 .* (2.0 * Spm .* (3.0 + 2.0*cos(x) + cos(2*x)) + Sop.* (2.0 + cos(x)));
    
    Sae(isnan(Sae)) = 1e30;
end