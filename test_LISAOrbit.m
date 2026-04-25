%
%   animation of LISAOrbit
%
% LDC params
theta_s                 =0.498944489240392;
phi_s                   =2.232796975920000;

Lt = 8.339102379953800437e+00;
shiftL = 3*Lt;
dur                     =62914560;
dt                      =15.0;
tVecl = 0:dt:dur;
dt_orbit = 61440;
tVecs_lisaorbit = 0:dt_orbit:dur;

[xt,yt,zt,n1s,n2s,n3s,R0,R1,R2,R3,q1,q2,q3,Lt] = func_LISAOrbit(tVecs_lisaorbit);

pic_num = 1;

for j=1:5:length(tVecs_lisaorbit)   %%%    for循环

    x = [xt(1,j),xt(2,j),xt(3,j),xt(1,j)];
    y = [yt(1,j),yt(2,j),yt(3,j),yt(1,j)];
    z = [zt(1,j),zt(2,j),zt(3,j),zt(1,j)];
    
    plot3(x,y,z,'b');
    hold on
    plot3(0,0,0,'mo','MarkerSize',30);
    hold off
    xlim([-500,500])
    ylim([-500,500])
    zlim([-4,4])
    grid on
    %view(j,18);     %%%%    移动视角
    pause(0.09);     %%%%    暂停时间
    
    F=getframe(gcf);
    I=frame2im(F);
    [I,map]=rgb2ind(I,256);
    if pic_num == 1
    imwrite(I,map,'./Figures/test_2.gif','gif','Loopcount',inf,'DelayTime',0.2);
    else
    imwrite(I,map,'./Figures/test_2.gif','gif','WriteMode','append','DelayTime',0.2);
    end
    pic_num = pic_num + 1;
end