
dt = 15;
T = 86400*365;
df = 1/T; 
fs = 1/dt;
fVec = 0:df:fs/2;

%model = 'mldc';
%model = 'Proposal';
model ='SciRDv1';
unit = 'relativeFrequency';
%unit = 'displacement';

Sx = noisepsd_X(fVec, model, unit);
Sae = noisepsd_AE(fVec,model,unit);
St = noisepsd_T(fVec,model,unit);


figure
plot(log10(fVec),log10(Sx),log10(fVec),log10(Sae),log10(fVec),log10(St));
legend('psd x','psd ae','psd t')


[Sx,Sae,St] = noisepsd_1(fVec,model,unit);

figure
plot(log10(fVec),log10(Sx),log10(fVec),log10(Sae),log10(fVec),log10(St));
legend('psd x','psd ae','psd t')