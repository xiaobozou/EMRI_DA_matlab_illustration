function freqVec = fftfreqvec1(nSmpls,fs)
%Generates the frequency vector for an FFT 
%F = FFTFREQVEC(N)
%N is the number of data samples. F is the set of positive and negative
%frequencies arranged in the same way as the FFT of any data vector
%containing N samples.

%Soumya D. Mohanty, Dec 2019

%dataLenInv = 1/(nSmpls-1);
%posFreq = (0:floor(nSmpls/2))*dataLenInv;
%negFreq = (1:(nSmpls-length(posFreq)))*dataLenInv;

df = fs/nSmpls;
%kNyq = floor(nSmpls/2) + 1;
%posFreq = (0:(kNyq-1))*df;
if mod(nSmpls,2)
    kNyq = floor(nSmpls/2);
    posFreq = (1:(kNyq))*df;
    negFreq = -flip(posFreq);
    freqVec = [0, posFreq, negFreq];
else
    kNyq = nSmpls/2-1;
    posFreq = (1:(kNyq))*df;
    negFreq = -flip(posFreq);
    freqVec = [0, posFreq, 0, negFreq];
end


