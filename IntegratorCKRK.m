function yend=IntegratorCKRK(ystart,dydt,tstart,dt,u,M,lam,aa)

a=zeros(1,7);
bb = zeros(7,6);
cc = zeros(1,7);
dc = zeros(1,7);


  a(3)=0.2; a(4)=0.3; a(5)=0.6; a(6)=1.0; a(7)=0.875;
  bb(3,2)=0.2;
  bb(4,2)= 3.0/40.0; bb(4,3)= 9.0/40.0;
  bb(5,2)=0.3; bb(5,3) = -0.9; bb(5,4) = 1.2;
  bb(6,2)=-11.0/54.0; bb(6,3)=2.5; bb(6,4)=-70.0/27.0; bb(6,5)=35.0/27.0;
  bb(7,2)=1631.0/55296.0; bb(7,3)=175.0/512.0; bb(7,4)=575.0/13824.0;
  bb(7,5)=44275.0/110592.0; bb(7,6)=253.0/4096.0;
  cc(2)=37.0/378.0; cc(4)=250.0/621.0; cc(5)=125.0/594.0; cc(7)=512.0/1771.0;
  dc(6)=-277.0/14336.0; dc(2)=cc(2)-2825.0/27648.0; dc(4)=cc(4)-18575.0/48384.0;
  dc(5)=cc(5)-13525.0/55296.0; dc(7)=cc(7)-0.25;

 ytemp = zeros(1,5); 
 ak2 = zeros(1,5); 
 ak3= zeros(1,5); 
 ak4 = zeros(1,5); 
 ak5 = zeros(1,5); 
 ak6 = zeros(1,5); 
 yend = zeros(1,5); 
 for i=1:5
      ytemp(i) = ystart(i) + bb(3,2)*dt*dydt(i);
 end
   %fprintf("RK step 1 \n");
 ak2 = Orbit_ODE(tstart+a(3)*dt,ytemp,u,M,lam,aa);
 
 for i=1:5
      ytemp(i) = ystart(i) + ( bb(4,2)*dydt(i) + bb(4,3)*ak2(i)  )*dt;
 end
    %fprintf("RK step 2 \n");
 ak3 = Orbit_ODE(tstart+a(4)*dt,ytemp,u,M,lam,aa);
 
  for i=1:5
      ytemp(i) = ystart(i) + ( bb(5,2)*dydt(i) + bb(5,3)*ak2(i) + bb(5,4)*ak3(i)  )*dt;
  end
     %fprintf("RK step 3 \n");
  ak4 = Orbit_ODE(tstart+a(5)*dt,ytemp,u,M,lam,aa);
  
  for i=1:5
      ytemp(i) = ystart(i) + ( bb(6,2)*dydt(i) + bb(6,3)*ak2(i) + bb(6,4)*ak3(i) + bb(6,5)*ak4(i) )*dt;
  end
     %fprintf("RK step 4 \n");
  ak5 = Orbit_ODE(tstart+a(6)*dt,ytemp,u,M,lam,aa);
  
  for i=1:5
      ytemp(i) = ystart(i) + ( bb(7,2)*dydt(i) + bb(7,3)*ak2(i) + bb(7,4)*ak3(i) + bb(7,5)*ak4(i) + bb(7,6)*ak5(i))*dt;
  end
     %fprintf("RK step 5 \n");
  ak6 = Orbit_ODE(tstart+a(7)*dt,ytemp,u,M,lam,aa);
  
  
  for i=1:5
      yend(i) = ystart(i) + ( cc(2)*dydt(i) + cc(4)*ak3(i)  + cc(5)*ak4(i) + cc(7)*ak6(i))*dt;
  end
  
  %fprintf("\n");
%fprintf("yend(1) = %.18e\n", yend(1));
%fprintf("yend(2) = %.18e\n", yend(2));
%fprintf("yend(3) = %.18e\n", yend(3));
%fprintf("yend(4) = %.18e\n", yend(4));
%fprintf("yend(5) = %.18e\n", yend(5));
%fprintf("\n");
end