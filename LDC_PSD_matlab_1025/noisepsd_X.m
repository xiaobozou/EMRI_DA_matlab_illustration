function Sx = noisepsd_X(fvec, model,unit)
   
clight = 299792458.0;    %# Speed of light in vacuum [m.s^-1] (CODATA 2014)
c = clight;
lisaL = 2.5e9;% # LISA's arm meters
lisaLT = lisaL/clight;% # LISA's armn in sec


    x = 2.0 * pi * lisaLT * fvec;

    [Spm, Sop] = lisanoises(fvec, model,unit);

    Sx = 16.0 * sin(x).^2 .* (2.0 * (1.0 + cos(x).^2) .* Spm + Sop);
    
    Sx(isnan(Sx)) = 1e30;

end