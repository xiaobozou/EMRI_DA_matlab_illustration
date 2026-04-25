function [Fp,Fc,k] = func_FpFc(eclipticLat,eclipticLong,t,n1,n2,n3)
% we are using HAK without psi
% beta = pi/2 - eclipticLat
% lambda = eclipticLong

%GW propagation
% LDC manual (10)
k = [-cos(eclipticLat)*cos(eclipticLong), -cos(eclipticLat)*sin(eclipticLong), -sin(eclipticLat)];
u = [sin(eclipticLong), -cos(eclipticLong), 0];
v = [-sin(eclipticLat)*cos(eclipticLong), -sin(eclipticLat)*sin(eclipticLong), cos(eclipticLat)];

% polarization tensor
PT_p = transpose(u)*u - transpose(v)*v; 
PT_c = transpose(u)*v + transpose(v)*u; 


% arm vector
%[xt,yt,zt,n1,n2,n3,R1,R2,R3,Lt] = func_LISAOrbit(t);
% n1~[3,N]

% detector response
Fp = zeros(3,length(t));
Fc = zeros(3,length(t));
% 
% for i=1:length(t)
%     
%     % detector tensor of single arm for moment i
%     D1 = n1(:,i)*transpose(n1(:,i));
%     D2 = n2(:,i)*transpose(n2(:,i));
%     D3 = n3(:,i)*transpose(n3(:,i));
%     
% 
% 
%     Fp(1,i) = sum(sum(D1.*PT_p)); 
%     Fp(2,i) = sum(sum(D2.*PT_p)); 
%     Fp(3,i) = sum(sum(D3.*PT_p));

%     Fc(1,i) = sum(sum(D1.*PT_c)); 
%     Fc(2,i) = sum(sum(D2.*PT_c)); 
%     Fc(3,i) = sum(sum(D3.*PT_c));
% end


for i=1:length(t)
    
    % detector tensor of single arm for moment i
%     D1 = n1(:,i)*transpose(n1(:,i));
%     D2 = n2(:,i)*transpose(n2(:,i));
%     D3 = n3(:,i)*transpose(n3(:,i));
    

    n1_t = transpose(n1(:,i));
    n2_t = transpose(n2(:,i));
    n3_t = transpose(n3(:,i));
        
    Fp(1,i) = n1_t*PT_p*n1(:,i); 
    Fp(2,i) = n2_t*PT_p*n2(:,i); 
    Fp(3,i) = n3_t*PT_p*n3(:,i); 
    
    
    Fc(1,i) = n1_t*PT_c*n1(:,i); 
    Fc(2,i) = n2_t*PT_c*n2(:,i); 
    Fc(3,i) = n3_t*PT_c*n3(:,i); 
end

end
