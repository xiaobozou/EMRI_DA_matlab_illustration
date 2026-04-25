function [Sx,Sae,St] = noisepsd_1(fVec,model,unit)

Sx = noisepsd_X(fVec, model,unit);
Sae = noisepsd_AE(fVec,model,unit);
St = noisepsd_T(fVec, model,unit);

end