function [kn,kR] = interp_knkR(k,t1,R1s,R2s,R3s,n1s,n2s,n3s,t2)

% 2022/03/28 interp knkR
N = length(t2);
kn = zeros(3,N);
kR = zeros(3,N);

kn1s = k*n1s;
kn2s = k*n2s;
kn3s = k*n3s;

F = griddedInterpolant(t1, kn1s);
kn(1,:) = F(t2);
F = griddedInterpolant(t1, kn2s);
kn(2,:) = F(t2);
F = griddedInterpolant(t1, kn3s);
kn(3,:) = F(t2);


kR1s = k*R1s;
kR2s = k*R2s;
kR3s = k*R3s;

F = griddedInterpolant(t1, kR1s);
kR(1,:) = F(t2);
F = griddedInterpolant(t1, kR2s);
kR(2,:) = F(t2);
F = griddedInterpolant(t1, kR3s);
kR(3,:) = F(t2);
end