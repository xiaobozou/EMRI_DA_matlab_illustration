function [hpSSB,hcSSB]=func_get_HAK_nm(theta_s,phi_s,theta_k,phi_k,ut,Mt,aa,lam,Dt,phit,nut,gamt,et,alpt,n,m)
% 2022/0120 change the input
% use ode solution as input, rather than params

% 2022/03/28
N = length(phit); % 2^22+1

% in Ap,As, using polar angles of thetaS
theta_s   = pi/2 - theta_s;

%% 　ComputeHarmonics
% tic;
Amp = ut/Dt;

cS= cos(theta_s);
sS= sin(theta_s);
cK= cos(theta_k);
sK= sin(theta_k);
cSp= cos(phi_s);
sSp= sin(phi_s);
cKp= cos(phi_k);
sKp= sin(phi_k);
clam=cos(lam);
slam=sin(lam);


Sn = cS*cK + sS*sK*cos(phi_s  - phi_k);
cX = Sn;

sX = sqrt( sS*sS*cK*cK - 2.0*sS*sSp*cK*cS*sK*sKp + ...
                cS*cS*sK*sK - 2.0*cS*sK*cKp*sS*cSp*cK + ...
                sS*sS*cSp*cSp*sK*sK*sKp*sKp - 2.0*sS*sS*cSp*sK*sK*sKp*sSp*cKp +...
                sS*sS*sSp*sSp*sK*sK*cKp*cKp);

Apc1 = ( -cK*cKp*sS*cSp - cK*sKp*sS*sSp + sK*cS )/(sX);
Aps1 = ( sKp*sS*cSp - cKp*sS*sSp )/(sX);
Apcn1 = ( cK*cS + sK*cKp*sS*cSp + sK*sKp*sS*sSp - cX)*clam/(sX*slam);
Apc2 = (sS*cSp*sKp - sS*sSp*cKp )*clam/(sX);
Aps2 = ( cK*cKp*sS*cSp + cK*sKp*sS*sSp - cS*sK )*clam/(sX);
Apcn2 = 0.0;


Aqc1 = ( sS*cSp*sKp - sS*sSp*cKp  )*cX/(sX);
Aqs1 = ( cK*cKp*sS*cSp + cK*sKp*sS*sSp - cS*sK )*cX/(sX);
Aqcn1 = 0.0;


Aqc2 = cX*clam*( cK*cKp*sS*cSp + cK*sKp*sS*sSp - sK*cS)/(sX);
Aqs2 = -cX*clam*sS*( sKp*cSp - cKp*sSp )/(sX);
Aqcn2 = -( cX*clam*clam*( cK*cS + sK*cKp*sS*cSp + sK*sKp*sS*sSp ) +...
                      1.- cX*cX - clam*clam )/(sX*slam);

Bp1c1 = 2.0*(Apc1*Apcn1 - Aqc1*Aqcn1 + Aqc2*Aqcn2 - Apc2*Apcn2);
Bp1c2 =  0.5*(Aps2*Aps2 - Aqc1*Aqc1  + Apc1*Apc1  - Aps1*Aps1 + ...
                          Aqc2*Aqc2 + Aqs1*Aqs1 - Apc2*Apc2 - Aqs2*Aqs2);
Bp1s1 = 2.0*(Aqs2*Aqcn2 - Aps2*Apcn2 - Aqs1*Aqcn1 + Aps1*Apcn1);
Bp1s2 = (Apc1*Aps1 + Aqc2*Aqs2 - Apc2*Aps2 - Aqc1*Aqs1);
Bp1cn = 0.5*(Apc1*Apc1 + Aps1*Aps1 - Aqc1*Aqc1 - Aqs1*Aqs1 - Apc2*Apc2 ...
                          + Aqc2*Aqc2 + Aqs2*Aqs2 - Aps2*Aps2) + Aqcn2*Aqcn2 - Aqcn1*Aqcn1 ...
                          + Apcn1*Apcn1 - Apcn2*Apcn2;

Bp2c1 = (Apcn1*Apc2 + Apc1*Apcn2 - Aqcn1*Aqc2 - Aqc1*Aqcn2);
Bp2c2 = 0.5*(Aqs1*Aqs2 - Aps1*Aps2 + Apc1*Apc2 - Aqc1*Aqc2);
Bp2s1 = (Aps1*Apcn2 + Apcn1*Aps2 - Aqcn1*Aqs2 - Aqs1*Aqcn2);
Bp2s2 = 0.5*( Apc1*Aps2 - Aqc1*Aqs2 + Aps1*Apc2 - Aqs1*Aqc2);
Bp2cn = 0.5*(Aps1*Aps2 - Aqs1*Aqs2 - Aqc1*Aqc2 + Apc1*Apc2) -Aqcn1*Aqcn2 + Apcn1*Apcn2;

Bc1c1 = (-Apc2*Aqcn2 - Apcn2*Aqc2 + Apc1*Aqcn1 + Apcn1*Aqc1);
Bc1c2 = 0.5*( Apc1*Aqc1 - Aps1*Aqs1 - Apc2*Aqc2 + Aps2*Aqs2);
Bc1s1 = (Apcn1*Aqs1 - Aps2*Aqcn2 + Aps1*Aqcn1 - Apcn2*Aqs2);
Bc1s2 = 0.5*(-Apc2*Aqs2 + Apc1*Aqs1 - Aps2*Aqc2 + Aps1*Aqc1);
Bc1cn = -Apcn2*Aqcn2 + Apcn1*Aqcn1 + 0.5*(Apc1*Aqc1 - Aps2*Aqs2 + Aps1*Aqs1 - Apc2*Aqc2);

Bc2c1 = (Aqc1*Apcn2 + Aqcn1*Apc2 + Apc1*Aqcn2 + Apcn1*Aqc2);
Bc2c2 = 0.5*( Apc1*Aqc2 - Aps1*Aqs2 + Aqc1*Apc2 - Aqs1*Aps2);
Bc2s1 = (Apcn1*Aqs2 + Aqs1*Apcn2 + Aps1*Aqcn2 + Aqcn1*Aps2);
Bc2s2 = 0.5*(Aqc1*Aps2 + Apc1*Aqs2 + Aqs1*Apc2 + Aps1*Aqc2);
Bc2cn = Aqcn1*Apcn2 + Apcn1*Aqcn2 + 0.5*(Apc1*Aqc2 + Aqs1*Aps2 +Aps1*Aqs2 + Aqc1*Apc2);

%------------------------
Ap = zeros(1,5);
Ac = zeros(1,5);
Ap(1) = 0.5*((Bp1c2+Bp2s2) - 1.0j*(Bp2c2-Bp1s2));  %# Aplus = Re( Ap*exp(1j*phi))
Ac(1) = 0.5*((Bc1c2+Bc2s2) - 1.0j*(Bc2c2-Bc1s2)); %# Across = Re( Ac*exp(1j*phi))
Ap(2) = 0.5*((Bp1c1+Bp2s1) - 1.0j*(Bp2c1-Bp1s1));
Ac(2) = 0.5*((Bc1c1+Bc2s1) - 1.0j*(Bc2c1-Bc1s1));
Ap(3) = Bp1cn - 1.0j*Bp2cn;
Ac(3) = Bc1cn - 1.0j*Bc2cn;
Ap(4) = 0.5*((Bp1c1-Bp2s1) - 1.0j*(Bp2c1+Bp1s1));
Ac(4) = 0.5*((Bc1c1-Bc2s1) - 1.0j*(Bc2c1+Bc1s1));
Ap(5) = 0.5*((Bp1c2-Bp2s2) - 1.0j*(Bp2c2+Bp1s2));
Ac(5) = 0.5*((Bc1c2-Bc2s2) - 1.0j*(Bc2c2+Bc1s2));

%----------------------- constant psi ---------------
up = cS*sK*cos(phi_s-phi_k) - cK*sS;
dw = sK*sin(phi_s-phi_k);
psi = atan2(up, dw);
fprintf("psi in harmonics comp =%f\n", psi);
cps2 = cos(2.0*psi);
sps2 = sin(2.0*psi);

% time = toc;
% fprintf("Ap,Ac  take %f sec\n",time);

%% single harmonics 
%---------------------Aps,Acs,Phi-----------------------
tic;
%phases = zeros(5, 5, sz_t);
%Aps    = zeros(5, 5, sz_t);
%Acs    = zeros(5, 5, sz_t);
% 2022/02/25 not store amplitude and phase of each harmonics, just use them
phases = zeros(1, N);
Aps    = zeros(1, N);
Acs    = zeros(1, N);

% time-dependent amplitude
Ampl = (2.0*pi*nut*Mt).^(2.0/3.0)*Amp;
%     e2 = et.*et;
%     e3 = e2.*et;
%     e4 = e2.*e2;
%     e5 = e3.*e2;
%     e6 = e4.*e2;
%     e7 = e5.*e2;

% in func_get_myHAK_3.m
et2 = et.*et;
et3 = et2.*et;
et4 = et2.*et2;
et5 = et3.*et2;
et6 = et4.*et2;
et7 = et5.*et2;

% time dependent ampltude
% Xp = zeros(6,N);
% % in func_get_myHAK_3.m
% Xp(1,:) = 0;
% Xp(2,:) = -3.*et + 13./8.*et3;
% Xp(3,:) = 2. - 5.*et2 + 23./8.*et4;
% Xp(4,:) = 3.*et - 57./8.*et3 + 321./64.*et5;
% Xp(5,:) = 4.*et2 - 10.*et4 + 101./12.*et6;
% Xp(6,:) = 125./24.*et3 - 5375./384.*et5 + 42125./3072.*et7;

nXp = zeros(6,N);
% in func_get_myHAK_3.m
nXp(1,:) = zeros(1,N);
nXp(2,:) = 1.0*(-3.*et + 13./8.*et3);
nXp(3,:) = 2.0*(2. - 5.*et2 + 23./8.*et4);
nXp(4,:) = 3.0*(3.*et - 57./8.*et3 + 321./64.*et5);
nXp(5,:) = 4.0*(4.*et2 - 10.*et4 + 101./12.*et6);
nXp(6,:) = 5.0*(125./24.*et3 - 5375./384.*et5 + 42125./3072.*et7);


gamt2 = 2.0*gamt;

hpS = zeros(1,N);
hcS = zeros(1,N);

hpSSB = zeros(1,N);
hcSSB = zeros(1,N);
%#### saving harmonics separately
%for n=1:5
%    %  error  nXp = n*Xp(n)
%    % 0105 bug
%    % nXp = n*Xp(n+1);  // bug
%    nXp = n*Xp(n+1,:);   //corrected 
%    for ai=1:5
%        phases(n, ai, :) =  phit*n + gma2 + (ai-3.0)*alpt; %### n, 2, ai : phi, gamma, alpha
%        Aps(n, ai, :)  = -0.5*Ampl.*nXp*Ap(ai);
%        Acs(n, ai, :)  = -Ampl.*nXp*Ac(ai);
%        
%        hpS = hpS + reshape(real(Aps(n, ai, :).*exp(1.0j*phases(n, ai, :))), 1,sz_t);
%        hcS = hcS + reshape(real(Acs(n, ai, :).*exp(1.0j*phases(n, ai, :))), 1,sz_t);
%    end
%end
        
        %nXp = n*Xp(n+1,:);
        
        phases =  n*phit + gamt2 + (m-3.0)*alpt;
        % old
        %Aps = -0.5*Ampl.*nXp*Ap(m);
        %Acs = -Ampl.*nXp*Ac(m);
        % duplicate calculation of Ampl.*nXp
        
        %hpS = hpS + real(Aps.*exp(1.0j*phases));
        %hcS = hcS + real(Acs.*exp(1.0j*phases));
        
        % 2022/03/28 update HAK
        Amp_nt = Ampl.*nXp(n+1,:).*exp(1.0j*phases);
        hpS = -0.5*real(Amp_nt*Ap(m));
        hcS = -1.0*real(Amp_nt*Ac(m));
       

% from source to SSB
tic;
hpSSB  =  hpS*cps2 + hcS*sps2;
hcSSB  = -hpS*sps2 + hcS*cps2;
time = toc;
fprintf("HAK from source to SSB take %f sec\n", time);

clear hpS;
clear hcS;
end
