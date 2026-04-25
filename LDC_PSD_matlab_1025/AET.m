function [A,E,T] = AET(X,Y,Z)

    A = (Z - X)/sqrt(2.0); 
    E = (X - 2.0*Y + Z)/sqrt(6.0);
    T = (X + Y + Z)/sqrt(3.0);
end