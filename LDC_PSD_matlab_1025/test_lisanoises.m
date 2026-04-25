
dt = 15;
T = 86400*365;
df = 1/T; 
fs = 1/dt;
fvec = 0:df:fs/2;

%model = 'mldc';
model = 'Proposal';
%unit = 'relativeFrequency';
unit = 'displacement';
[Spm, Sop]=lisanoises(fvec,model,unit);