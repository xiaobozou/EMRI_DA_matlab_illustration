function LLR_TD = fitness_10D(x,params)

    ut       =     params.ut;
    Mt       =     params.Mt;
    lam     =     params.lam;
    aa      =     params.aa;
    e0      =     params.e0;
    nu0     =     params.nu0;
    theta_s =     params.theta_s;
    phi_s   =     params.phi_s;
    theta_k =     params.theta_k;
    phi_k   =     params.phi_k;
    Dt      =     params.Dt;

    phi0 = x(1);
    gam0 = x(2);
    alp0 = x(3);


    tVecl =     params.tVecl;
    tVecs_ode =     params.tVecs_ode;
    tVecs_lisaorbit =     params.tVecs_lisaorbit;

    Sae2 =      params.Sae2;
    Adw2 = params.Adw2;
    Edw2 = params.Edw2;

    dur = params.dur;
    dt = params.dt;
    dt_ode = params.dt_ode;
    id_cut = params.id_cut;

    Nl = length(tVecl);
    N_fft = dur/dt;    
    fs = 1.0/dt;
    df = fs/N_fft;
    %%  fidutial orbit and response 
    [xts,yts,zts,n1s,n2s,n3s,R0s,R1s,R2s,R3s,q1s,q2s,q3s,Lt]  = func_LISAOrbit(tVecs_lisaorbit);
    [Fps,Fcs,k] = func_FpFc(theta_s,phi_s,tVecs_lisaorbit,n1s,n2s,n3s);
    
    %% interp response 
    [Fp,Fc] = interp_FpFc(tVecs_lisaorbit,Fps,Fcs,tVecl);
    [kn,kR] = interp_knkR(k,tVecs_lisaorbit,R1s,R2s,R3s,n1s,n2s,n3s,tVecl);
    
    % -----------------  ODE->HAK->TDI -----------------%
    do_shift = 1; % shift template XYZ, start from 42
    %do_shift = 0; 
    %% orbital ode
    %  check_ode_start
    start_ok = check_ode_start(e0,nu0,Mt);
    fprintf("start_ok = %d\n",start_ok);

    if start_ok == 0
        fprintf("check_ode_start   fail \n");
    end
        
    % solve ode  in dt_ode
    [phit1,nut1,gamt1,et1,alpt1,len] = solve_ode(phi0,nu0,gam0,e0,alp0,tVecs_ode,ut,Mt,lam,aa);
    fprintf("for dt_ode, len = %d,  tVecs_ode(len) = %f\n", len, tVecs_ode(len));


    % interp ode
    sz_t = (len-1)*dt_ode/dt+1;
    fprintf("for dt, sz_t = %d,  tVecl(sz_t) = %f\n",  sz_t, tVecl(sz_t));
    
    phit = zeros(1,Nl);
    nut  = zeros(1,Nl);
    gamt = zeros(1,Nl);
    et   = zeros(1,Nl);
    alpt = zeros(1,Nl);
    [phit,nut,gamt,et,alpt] = interp_ode(tVecs_ode, phit1,nut1,gamt1,et1,alpt1,tVecl);
   
    %% HAK
    % take interolated ode tp 14D hphc generation
    % 25 harmonics
    num_harmonics = 25;
    % 10 harmonics
    %num_harmonics = 10;
    % 5 harmonics
    %num_harmonics = 5; 
    % 3 harmonics
    %num_harmonics = 3;  
    [hp,hc]=func_get_HAK(theta_s,phi_s,theta_k,phi_k,ut,Mt,aa,lam,Dt,phit,nut,gamt,et,alpt,num_harmonics);

    %% TDI 1.0
    slr_X = [1,-3, 2; 2, 3, 1;1,2, 3;3,-2, 1;1,2, 3;3,-2, 1;1,-3, 2;2, 3, 1];
    slr_Y = [2,-1, 3; 3, 1, 2;2,3, 1;1,-3, 2;2,3, 1;1,-3, 2;2,-1, 3;3, 1, 2];
    slr_Z = [3,-2, 1; 1, 2, 3;3,1, 2;2,-1, 3;3,1, 2;2,-1, 3;3,-2, 1;1, 2, 3];
    slr_idx = [slr_X;slr_Y;slr_Z];
    shiftL = [3*Lt,2*Lt,1*Lt,0*Lt,3*Lt,2*Lt,1*Lt,0*Lt];
    shiftL1 = [shiftL,shiftL,shiftL];

    tic;
    [X,Y,Z] = func_XYZ_Mohanty(hp,hc,Fp,Fc,kn,kR,tVecl,Lt,slr_idx,shiftL1);
    time = toc;
    fprintf("total XYZ cost = %f sec\n",time);

    if do_shift == 1
        X = circshift(X,-41);
        Y = circshift(Y,-41);
        Z = circshift(Z,-41);
    end

   [A,E,~] = AET(X(1:N_fft),Y(1:N_fft),Z(1:N_fft));
    Af  = fft(A)*dt;   Ef = fft(E)*dt; 
    innerprod_hh_A_FD = innerproduct_fre_2(Af,Af,Sae2,df,id_cut);
    innerprod_hh_E_FD = innerproduct_fre_2(Ef,Ef,Sae2,df,id_cut); 
    
    innerprod_dh_A_TD = 2*sum(A.*Adw2)*dt;
    innerprod_dh_E_TD = 2*sum(E.*Edw2)*dt;
    LLR2_TD =  (innerprod_dh_A_TD+innerprod_dh_A_TD)^2/(innerprod_hh_A_FD+innerprod_hh_A_FD);
    LLR_TD = sqrt(LLR2_TD);
    fprintf("do_shift=%d, LLR TD =%f\n",do_shift, LLR_TD);

end
