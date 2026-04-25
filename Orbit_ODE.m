function dydt=Orbit_ODE(t,y,u,M,lam,aa)

% 
twopi = 2*pi;
uM3 = u/(M*M*M);
uM2 = u/(M*M);
p1=  twopi*M*y(2);

e = y(4);
e2 =  e*e;
 
e4 = e2*e2;
e6 = e4*e2;
p2 = 1-e2;

ip2 = 1/p2 ;
ip2_2 = ip2*ip2;
s_ip2 = sqrt(ip2);

clam = cos(lam);
p3 = twopi*y(2);

% ###################################################### 
% in matlab, pow(fr, 2./3.) =  8.759259828850035623e-02
% in c,      pow(fr, 2./3.) =  8.759259828850034235e-02

% pow(fr, 8./3.) equals.

% in matlab, pow(fr, 11./3.) = 1.526056982544237838e-06
% in c,      pow(fr, 11./3.) = 1.526056982544238050e-06

% pow(fr, 1./3.)*pow(fr, 1./3.) equals.
% pow(fr, 4.)/pow(fr, 1./3.) equals.
% pow(fr, 3.)/pow(fr, 1./3.) equals.

% that result in the ode difference.
% use pow(fr, 1./3.)*pow(fr, 1./3.) to replace pow(fr, 2./3.)
% use pow(fr, 3.)/pow(fr, 1./3.) to replace    pow(fr, 8./3.)
% use pow(fr, 4.)/pow(fr, 1./3.) to replace    pow(fr, 11./3.)
%

p4 = p1^(2/3);          %p4 = p1^(1/3)*p1^(1/3);
p5 = p1^(8/3);          %p5 = p1^(3)/p1^(1/3); 
p6 = p1^(11/3);         %p6 = p1^(4)/p1^(1/3);
% ###################################################### %


dydt = zeros(5,1);

% y(1) = Phi,
dydt(1) = twopi*y(2);

% y(2) =  nu,
dydt(2) =9.6/pi*uM3*p6*ip2_2*ip2_2*s_ip2*...
    (  (1+73/24*e2+37/96*e4)*p2 + p4*(1273/336-2561/224*e2-3885/128*e4-13147/5376*e6) ...
    -p1*aa*clam*s_ip2*(73/12 +1211/24*e2+3143/96*e4 +65/64*e6) );  

% y(3) =  gam, 
% s_ip2*ip2 -> ip2*s_ip2
dydt(3) =3*p3*p4*ip2*(1+0.25*p4*ip2*(26-15*e2)) -6*p3*clam*aa*p1*ip2*s_ip2;

% y(4) = e,  8*16705->133640,  12*9082-> 108984 
dydt(4) = -e/15*uM2*s_ip2*ip2*ip2_2*p5*...
    ((304+121*e2)*p2*(1+12*p4)-1/56*p4*(133640+108984*e2-25211*e4) ) ...
    + e*uM2*aa*clam*p6*ip2_2*ip2_2*(1364/5 +5032/15*e2+26.3*e4);

% y(5) =  alp, ip2*s_ip2 -> s_ip2*ip2
dydt(5) = 2*p3*aa*p1*s_ip2*ip2;

%fprintf("\n");
%fprintf("dydt(1) = %.18e\n", dydt(1));
%fprintf("dydt(2) = %.18e\n", dydt(2));
%fprintf("dydt(3) = %.18e\n", dydt(3));
%fprintf("dydt(4) = %.18e\n", dydt(4));
%fprintf("dydt(5) = %.18e\n", dydt(5));
%fprintf("\n");
end


