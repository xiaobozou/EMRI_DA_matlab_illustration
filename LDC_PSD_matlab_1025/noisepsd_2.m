function [Sx2,Sae2,St2] = noisepsd_2(Nt,fs,model,unit)

% model = 'mldc';
% %model = 'Proposal';
% unit = 'relativeFrequency';
df = fs/Nt;

% two-sided psd
if mod(Nt,2)
    kNyq = floor(Nt/2);
    % positive frequency array
    fVec = (1:kNyq)*df;
    % one-sided psd
    Sx = noisepsd_X(fVec, model, unit); 
    Sae = noisepsd_AE(fVec,model,unit); 
    St = noisepsd_T(fVec,model,unit);  
    % two-sided psd
    Sx2 = [1e30, Sx, flip(Sx)];
    Sae2 =[1e30, Sae, flip(Sae)];
    St2 = [1e30, St, flip(St)];    
else
    % positive frequency array
    kNyq = floor(Nt/2)-1;
    fVec = (1:kNyq)*df;
    % one-sided psd
    Sx = noisepsd_X(fVec, model, unit); 
    Sae = noisepsd_AE(fVec,model,unit); 
    St = noisepsd_T(fVec,model,unit);  
    % two-sided psd
    Sx2 = [1e30,  Sx,   1e30, flip(Sx)];
    Sae2 =[1e30,  Sae,  1e30, flip(Sae)];
    St2 = [1e30,  St,   1e30, flip(St)];  
end

end