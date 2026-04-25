T = 86400*365;
dt = 15;
Nt = round(T/dt);
fs = 1/dt;

model = 'mldc';
% %model = 'Proposal';
unit = 'relativeFrequency';
[Sx2,Sae2,St2] = noisepsd_2(Nt,fs,model,unit);
freqVec = fftfreqvec1(Nt,fs);

figure
plot(freqVec,log10(Sx2))
plot(log10(Sx2))