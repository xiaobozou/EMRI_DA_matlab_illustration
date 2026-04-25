function y=innerProd_fre_2SidedPSD(a,b,S,df)
    % two sided PSD
    %y = 2*real(sum(a.*conj(b)./S))*df;
    y = real(sum( (a.*conj(b)+ b.*conj(a))./S))*df;
end