function start_ok = check_ode_start(e0,nu0,Mt)
%% Two constraints to check whether given e0,nu0 can start ode sucessfully
    % 2022/01/20
    
    
    % 1st
    % BC04 equ(59) give Schwarzschild plunge condition, the start nu0 should
    % not reach and exceed plunge
    
    % nu_lso of given e0,
    nu_lso = ((1-e0*e0)/(6+2*e0))^(1.5)/(2*pi*Mt);
    
    
    % 2st
    % also from BC04 equ(59), inspired by Mohanty,the start e0 should
    % not reach and exceed plunge
    %nu_cut_off = 1e-3;
    nu_cut_off =  ((1/6)^(1.5))/(2*pi*Mt);
    Mp = (2*pi*Mt*nu_cut_off)^(2/3);
    e_threshold = sqrt(1-6*Mp+Mp^2) - Mp;



    if nu0 > nu_lso | e0 < e_threshold
        start_ok = 0;
    else
        start_ok = 1;
    end
    
end
