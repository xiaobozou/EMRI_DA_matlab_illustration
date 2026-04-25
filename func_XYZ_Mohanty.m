function [X,Y,Z] = func_XYZ_Mohanty(hp,hc,Fp,Fc,kn,kR,tVecl,Lt,slr_idx,shiftL1)

    % XYZ from Mohanty
    % 0405/2021
    % 1025/2021 modified
    % 0120/2022 modified
    
    Nl = length(tVecl);
        
    % XYZ
    %
    % --------------------------------------------------------------
    % 2022/03/28 cost 3.8 sec here, but 24 evluation of single y takes
    % 1.5 sec, there are 1.2 sec for self cost, now give up 24-D
    % y,using multile 1-D y instead
    % ----------------------------------------------------------------
    %
    
    X = zeros(1,Nl);
    Y = zeros(1,Nl);
    Z = zeros(1,Nl);

    id = 0;
    for i=1:24
        if i==1 | i==9 | i== 17
            temp = zeros(1,Nl);
            id = 1;
        end
        
        if id ==1| id==2 | id==3 | id==4
            y =   func_compute_yslr_2(hp,hc,Fp,Fc, kn,kR,tVecl,slr_idx(i,1),slr_idx(i,2),slr_idx(i,3), shiftL1(i),Lt);
            temp = temp + y;
            id = id+1;
        elseif id ==5| id==6 | id==7 | id==8
            y =   func_compute_yslr_2(hp,hc,Fp,Fc, kn,kR,tVecl,slr_idx(i,1),slr_idx(i,2),slr_idx(i,3), shiftL1(i),Lt);
            temp = temp - y;
            id = id+1;
        else
            fprintf("...error id ...\n");
        end    
    
        if i == 8
            X = temp;
        elseif i==16
            Y = temp;
        elseif i==24
            Z = temp;
        end
    end
       
    % X(1:sz_t) = sum(y(1:4,:),1) - sum(y((1:4)+4,:),1); 
    % Y(1:sz_t)  = sum(y((1:4)+8,:),1) - sum(y((1:4)+12,:),1); 
    % Z(1:sz_t)  = sum(y((1:4)+16,:),1) - sum(y((1:4)+20,:),1); 
end
