function y = func_compute_yslr_2(hp,hc,Fp,Fc,kn,kR,tVecl,s,l,r,shiftL,Lt)
% single arm response 
% written 0426
% xiaobo
% fft + interp


% single arm response 
% written 1225
% xiaobo
% no fft, all interp


% 2022/01/20 update
lp = abs(l);

temp = Fp(lp,:).*hp + Fc(lp,:).*hc; 
F = griddedInterpolant(tVecl,temp);

%% Phi1
Phi1 = F(tVecl -kR(s,:) -Lt -shiftL);

%% Phi2
Phi2 = F(tVecl -kR(r,:)-shiftL);

% 0816 I forgot the sign before kn, ZXB
if  l>0
    y = 0.5./(1-kn(l,:)).*(Phi1-Phi2);
elseif l<0
    y = 0.5./(1+kn(abs(l),:)).*(Phi1-Phi2);
end


% 2022/12/24 test spikes in the short end transient of XYZ
%fprintf("%d\n",length(tvec)); % pass
% fprintf("slr=%d%d%d\n",s,l,r);  % pass
% for i= (length(tvec)- 5): length(tvec)
%      fprintf("%d   %.18e   %.18e   %.18e\n",i,Phi1(i),Phi2(i),y(i));
% end
end
