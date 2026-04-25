function St = noisepsd_T(fvec,model,unit)


clight = 299792458.0;    %# Speed of light in vacuum [m.s^-1] (CODATA 2014)
c = clight;
lisaL = 2.5e9;% # LISA's arm meters
lisaLT = lisaL/clight;% # LISA's armn in sec

    x = 2.0 * pi * lisaLT * fvec;

    [Spm, Sop] = lisanoises(fvec, model,unit);

    St = 16.0 * Sop .* (1.0 - cos(x)) .* sin(x).^2 + 128.0 * Spm .* sin(x).^2 .* sin(0.5*x).^4;
    
    St(isnan(St)) = 1e30;
end
