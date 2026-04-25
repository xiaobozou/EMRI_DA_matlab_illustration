function [xt,yt,zt,n1,n2,n3,R0,R1,R2,R3,q1,q2,q3,Lt] = func_LISAOrbit(t)
%
% 2022/11/13 find R1,R2,R3 differ between matlab and c.
% which result from  constant differnet  from LDC and c.
%
N = length(t);
xt = zeros(3,N); yt = zeros(3,N); zt = zeros(3,N);
lambda = 0; kapa = 0;

% ####################################################################
% 2022/11/13 old constant, which result in XYZ/snr/loglike differences
% between c and matlab

% % year = 86400*365; 
% % c = 299792458;
% % omega = 2*pi/year;  % LDC_omega
% % alpha_t = omega*t + kapa;
% % a = 1.5e11; % meter
% % L = 2.5e9;  % meter
% % Lt = L/c; 
% % e = L/(2*sqrt(3)*a); 
% % a = a/c;
% ####################################################################
% 2022/11/13 matlab using LDC's constant same as Constant.cc and
% LISAConstant.py
% constant value used in c

%const double LDC_year = 31558149.763545600;
%const double LDC_omega =1.9923849908611068e-07;
%const double LDC_a = 149597870700.0; //% meter
%const double LDC_L = 2.5e9;
%const double LDC_Lt = 8.3391023799538;
%//const double LDC_e = LDC_L/(2.0*sqrt(3.0)*LDC_a); 
%const double LDC_e =0.004811252243246882;
%//const double LDC_at = LDC_a/LDC_c;
%const double LDC_at = 499.00478383615643;
year = 31558149.763545600;
c = 299792458;
omega =1.990986592768378515e-07;
alpha_t = omega*t + kapa;
%a = 149597870700.0; % in meter
L = 2.5e9;
Lt = 8.339102379953800437e+00;
e  = 4.824185218078991429e-03;
a  = 4.990047838361564345e+02; % in sec
% ####################################################################
% xt yt zt
for i =1:3
   beta = (i-1)*2*pi/3 + lambda; 
   xt(i,:) = a*cos(alpha_t) + a*e*(sin(alpha_t).*cos(alpha_t)*sin(beta) - (1+sin(alpha_t).^2 )*cos(beta));
   yt(i,:) = a*sin(alpha_t) + a*e*(sin(alpha_t).*cos(alpha_t)*cos(beta) - (1+cos(alpha_t).^2 )*sin(beta));
   zt(i,:) = -sqrt(3)*a*e*cos(alpha_t-beta);
end


% ###################################################################
% n1n2n3
%n1 = zeros(N,3); n2 = zeros(N,3); n3 = zeros(N,3);

% 2022/11/14 following LDC LW_simple.py
%     // in LW_simple.py, arm = 2.5e9, self.armt = LDC_Lt, constant.
%     //        n1 = np.array([x[1,:]-x[2,:], y[1,:]-y[2,:], z[1,:] - z[2,:]])/self.armt
%     //        n2 = np.array([x[2,:]-x[0,:], y[2,:]-y[0,:], z[2,:] - z[0,:]])/self.armt
%     //        n3 = np.array([x[0,:]-x[1,:], y[0,:]-y[1,:], z[0,:] - z[1,:]])/self.armt

% LDC n1 n2 n3
n1 = [(xt(2,:) -xt(3,:));  (yt(2,:) -yt(3,:));  (zt(2,:) -zt(3,:))]/Lt;
n2 = [(xt(3,:) -xt(1,:));  (yt(3,:) -yt(1,:));  (zt(3,:) -zt(1,:))]/Lt;
n3 = [(xt(1,:) -xt(2,:));  (yt(1,:) -yt(2,:));  (zt(1,:) -zt(2,:))]/Lt;

% 0815 correction
% as Mohanty
% armlen = sqrt( (xt(2,:) -xt(3,:)).^2  +  (yt(2,:) -yt(3,:)).^2 + (zt(2,:) -zt(3,:)).^2 );
% n1 = [(xt(2,:) -xt(3,:))./armlen;  (yt(2,:) -yt(3,:))./armlen;  (zt(2,:) -zt(3,:))./armlen];
% armlen = sqrt( (xt(3,:)  -xt(1,:)).^2 +  (yt(3,:) -yt(1,:)).^2 + (zt(3,:) -zt(1,:)).^2 );
% n2 = [(xt(3,:) -xt(1,:))./armlen;  (yt(3,:) -yt(1,:))./armlen;  (zt(3,:) -zt(1,:))./armlen];
% armlen = sqrt( (xt(1,:)  -xt(2,:)).^2 +  (yt(1,:) -yt(2,:)).^2 + (zt(1,:) -zt(2,:)).^2 );
% n3 = [(xt(1,:) -xt(2,:))./armlen;  (yt(1,:) -yt(2,:))./armlen;  (zt(1,:) -zt(2,:))./armlen];
% ###################################################################

% LDC R1 R2 R3
%R1 = [(xt(1,:)-a*cos(alpha_t));  (yt(1,:)-a*sin(alpha_t));  zt(1,:)];
%R2 = [(xt(2,:)-a*cos(alpha_t));  (yt(2,:)-a*sin(alpha_t));  zt(2,:)];
%R3 = [(xt(3,:)-a*cos(alpha_t));  (yt(3,:)-a*sin(alpha_t));  zt(3,:)];
R0 = [a*cos(alpha_t); a*sin(alpha_t); zeros(1,length(t))];

% 0815 correction
% as Mohanty
R1 = [xt(1,:);  yt(1,:);  zt(1,:)];
R2 = [xt(2,:);  yt(2,:);  zt(2,:)];
R3 = [xt(3,:);  yt(3,:);  zt(3,:)];
R0 = [a*cos(alpha_t);   a*sin(alpha_t);  zeros(1,length(t))];

q1 = [xt(1,:)-a*cos(alpha_t);   yt(1, :)-a*sin(alpha_t);     zt(1, :)];
q2 = [xt(2,:)-a*cos(alpha_t);   yt(2, :)-a*sin(alpha_t);     zt(2, :)];
q3 = [xt(3,:)-a*cos(alpha_t);   yt(3, :)-a*sin(alpha_t);     zt(3, :)];

end
