function y=innerproduct_fre_2(af2,bf2,Sn2,df,id_cut)

     % old
     %y = 2*real(sum(af2(id_cut:(end-id_cut)).*conj(bf2(id_cut:(end-id_cut)))./Sn2(id_cut:(end-id_cut))))*df;

     % take care of the half compex of fft(real A) when fitering LISAband
     % 2022/11/05
     % the length is N_fft = 2^22
     % new         
     af2(1) = 0;
     af2(2:id_cut) = 0;
     af2((end-(id_cut-2)):end) = 0;

     
     bf2(1) = 0;
     bf2(2:id_cut) = 0;
     bf2((end-(id_cut-2)):end) = 0;
     
     y = 2*real(sum(af2.*conj(bf2)./Sn2))*df; 
end
