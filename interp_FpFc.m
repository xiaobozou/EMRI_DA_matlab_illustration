function [Fp2,Fc2]=interp_FpFc(t1,Fp1,Fc1,t2)

N = length(t2);
Fp2 = zeros(3,N);
Fc2 = zeros(3,N);


F = griddedInterpolant(t1, Fp1(1,:));
Fp2(1,:) = F(t2);
F = griddedInterpolant(t1, Fp1(2,:));
Fp2(2,:) = F(t2);
F = griddedInterpolant(t1, Fp1(3,:));
Fp2(3,:) = F(t2);


F = griddedInterpolant(t1, Fc1(1,:));
Fc2(1,:) = F(t2);
F = griddedInterpolant(t1, Fc1(2,:));
Fc2(2,:) = F(t2);
F = griddedInterpolant(t1, Fc1(3,:));
Fc2(3,:) = F(t2);
end