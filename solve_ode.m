function [phit,nut,gamt,et,alpt,len] = solve_ode(phit0,nut0,gamt0,et0,alpt0,tVecs_ode,ut,Mt,lam,aa)

% the ode solutions vectors are always have length N = 2^22
% it will increase cost in ode interploation and waveform generation (longer vectors)
% but will decease more cost in fft (small prime factorization)
% 2022/02/25 update this function

deltat = tVecs_ode(2) - tVecs_ode(1);
N = length(tVecs_ode);

len=1; %  (len-1) is number of steps that ode evolve till plunge/ 2022/03/22

% 0105
%N = ceil(T/deltat);
%tVec = (0:N-1)*deltat; 
%tVec = zeros(1,N); 

phit = zeros(1,N);
nut  = zeros(1,N);
gamt = zeros(1,N);
et   = zeros(1,N);
alpt = zeros(1,N);

y = [phit0,nut0,gamt0,et0,alpt0];
phit(1) = y(1); 
nut(1) =  y(2); 
gamt(1) = y(3); 
et(1) =   y(4);  
alpt(1) = y(5); 

%  fprintf("------------------i = 0,  ti = 0---------------------------\n");
%fprintf("y(1) = %.18e\n", y(1));
%fprintf("y(2) = %.18e\n", y(2));
%fprintf("y(3) = %.18e\n", y(3));
%fprintf("y(4) = %.18e\n", y(4));
%fprintf("y(5) = %.18e\n", y(5));
%fprintf("\n");
% yi
%
%　2022/11/06 to solve XYZ spike, snr/loglike difference
%for i=1:N-1
for i=2:N
    %i;
    %fprintf("-----------------i = %d, ti = %.1f-------------------------- \n", i,ti);
    dydt = Orbit_ODE(tVecs_ode(i-1),y,ut,Mt,lam,aa);
    %vpa(y);
    %vpa(dydt);
    
    yend=IntegratorCKRK(y,dydt,tVecs_ode(i-1),deltat,ut,Mt,lam,aa);
    %vpa(yend);
    
    
    nu_lso = ((1-yend(4)*yend(4))/(6+2*yend(4)))^(1.5)/(2*pi*Mt);
    nu_i =  yend(2);
    %　2022/11/06 to solve XYZ spike, snr/loglike difference
    %if nu_i < nu_lso
    if nu_i <= nu_lso    
        %tVec(i+1) = ti + deltat;
        phit(i) = yend(1);
        nut(i) =  yend(2);
        gamt(i) = yend(3);
        et(i)  =  yend(4);
        alpt(i) = yend(5);
        % 0717
        y = yend;
        len = len+1;
    else
         break;    
    end
    
end

% throw off zeros segment after plunge
% tVec = tVec(1:len);
% phit = phit(1:len);
% nut = nut(1:len);
% gamt = gamt(1:len);
% et = et(1:len);
% alpt = alpt(1:len);

end

