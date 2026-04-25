%
% test FpFc,knkR
% 2026/04/23
% Xiaobo Zou

clc
clear

dur                     =62914560;
theta_s                 =0.498944489240392;
phi_s                   =2.232796975920000;
%%   dt = 15.0
dt                      =15.0;
Nt = round(dur/dt)+1;
tVecl_1 = (0:Nt-1)*dt;  % dense time label

[xtl,ytl,ztl,n1l,n2l,n3l,R1l,R2l,R3l,Lt]  = func_LISAOrbit(tVecl_1);
[Fp_a,Fc_a,k] = func_FpFc(theta_s,phi_s,tVecl_1,n1l,n2l,n3l);

kR1 = k*R1l;
kR2 = k*R2l;
kR3 = k*R3l;
kR_a = [kR1;kR2;kR3];

kn1 = k*n1l;
kn2 = k*n2l;
kn3 = k*n3l;
kn_a = [kn1;kn2;kn3];

%%   dt = 1536000 sec (~1 day) and interp
dt_orbit                      =61440;
Nt = round(dur/dt_orbit)+1;
tVecl_2 = (0:Nt-1)*dt_orbit;  % loose time label

[xts,yts,zts,n1s,n2s,n3s,R1s,R2s,R3s,Lt]  = func_LISAOrbit(tVecl_2);

[Fps,Fcs,k] = func_FpFc(theta_s,phi_s,tVecl_2,n1s,n2s,n3s);

% interp
[Fp_b,Fc_b]=interp_FpFc(tVecl_2,Fps,Fcs,tVecl_1);

% interp
[kn_b,kR_b] = interp_knkR(k,tVecl_2,R1s,R2s,R3s,n1s,n2s,n3s,tVecl_1);



%% plot
figure(1)
subplot(4,3,1)
plot(tVecl_1,Fp_a(1,:))
set(gca,'FontSize',20)  
xlabel('t (sec)')
ylabel('Fp')
legend('no interp')
subplot(4,3,2)
plot(tVecl_1,Fp_b(1,:))
set(gca,'FontSize',20)  
xlabel('t (sec)')
ylabel('Fp')
legend('interp')
subplot(4,3,3)
plot(tVecl_1,Fp_a(1,:)-Fp_b(1,:))
set(gca,'FontSize',20)  
xlabel('t (sec)')
ylabel('Fp')
legend('diff')

subplot(4,3,4)
plot(tVecl_1,Fc_a(1,:))
set(gca,'FontSize',20)  
xlabel('t (sec)')
ylabel('Fc')
legend('no interp')
subplot(4,3,5)
plot(tVecl_1,Fc_b(1,:))
set(gca,'FontSize',20)  
xlabel('t (sec)')
ylabel('Fc')
legend('interp')
subplot(4,3,6)
plot(tVecl_1,Fc_a(1,:)-Fc_b(1,:))
set(gca,'FontSize',20)  
xlabel('t (sec)')
ylabel('Fc')
legend('diff')

subplot(4,3,7)
plot(tVecl_1,kn_a(1,:))
set(gca,'FontSize',20)  
xlabel('t (sec)')
ylabel('kn')
legend('no interp')
subplot(4,3,8)
plot(tVecl_1,kn_b(1,:))
set(gca,'FontSize',20)  
xlabel('t (sec)')
ylabel('kn')
legend('interp')
subplot(4,3,9)
plot(tVecl_1,kn_a(1,:)-kn_b(1,:))
set(gca,'FontSize',20)  
xlabel('t (sec)')
ylabel('kn')
legend('diff')

subplot(4,3,10)
plot(tVecl_1,kR_a(1,:))
set(gca,'FontSize',20)  
xlabel('t (sec)')
ylabel('kR')
legend('no interp')
subplot(4,3,11)
plot(tVecl_1,kR_b(1,:))
set(gca,'FontSize',20)  
xlabel('t (sec)')
ylabel('kR')
legend('interp')
subplot(4,3,12)
plot(tVecl_1,kR_a(1,:)-kR_b(1,:))
set(gca,'FontSize',20)  
xlabel('t (sec)')
ylabel('kR')
legend('diff')

screenSize = get(0, 'ScreenSize'); % Get screen resolution [1,4](@ref)  
set(gcf, 'OuterPosition', screenSize); % Cover entire screen  
sgt = sgtitle('Test interpolation in detector response, 2026/04/23','Color','k','FontSize',30,'Interpreter', 'latex');
saveas(gcf,'./Figures/test_2.png')