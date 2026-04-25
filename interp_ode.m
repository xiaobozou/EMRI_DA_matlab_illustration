function [phit2,nut2,gamt2,et2,alpt2] = interp_ode(t1,phit1,nut1,gamt1,et1,alpt1,t2)
%% first solve ode in detlat = 15360, 
% if len == 1, the ode fails, 
% else interp the ode
% 2022/01/20

% 2022/02/25
%t2 = (0:2^22-1)*15;
%len2 = 2^22;

% interp phi
F = griddedInterpolant(t1, phit1);
phit2 = F(t2);


% interp nu
F = griddedInterpolant(t1,nut1);
nut2 = F(t2);


% interp gam
F = griddedInterpolant(t1,gamt1);
gamt2 = F(t2);


% interp e
F = griddedInterpolant(t1, et1);
et2 = F(t2);


% interp alp
F = griddedInterpolant(t1,alpt1);
alpt2 = F(t2);


%len2 = length(t2);
end