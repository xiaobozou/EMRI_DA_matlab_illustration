function fVec =fVec1(Nt,fs)

df = fs/Nt;
%fVec = 0:df:fs/2
    if mod(Nt/2)
        knyq = (Nt+1)/2;
        fVec=(0:knyq-1)*df;
    else
        knyq = Nt/2+1;
        fVec=(0:knyq-1)*df;
    end
    
%knyq = floor(Nt/2)+1;   
%fVec =(0:knyq-1)*df;     
end