function q = quaternion(C)
%This function was obtained from the book: Atmospheric and Space Flight
%Dynamics
%TEWARI, A. Atmospheric and Space Flight Dynamics. Boston, MA: Birkhauser,
%2007. ISBN: 978-0-8176-4437-6
T = trace(C);
qsq = [1+2*C(1,1)-T;1+2*C(2,2)-T;1+2*C(3,3)-T;1+T]/4;
[x,i] = max(qsq);
q = [0,0,0,0];
if i == 4
    q(4) = sqrt(x);
    q(1) = (C(2,3)-C(3,2))/(4*q(4));
    q(2) = (C(3,1)-C(1,3))/(4*q(4));
    q(3) = (C(1,2)-C(2,1))/(4*q(4));
end
if i == 3
    q(3) = sqrt(x);
    q(1) = (C(1,3)+C(3,1))/(4*q(3));
    q(2) = (C(3,2)+C(2,3))/(4*q(3));
    q(4) = (C(1,2)-C(2,1))/(4*q(3));
end
if i == 2
    q(2) = sqrt(x);
    q(1) = (C(1,2)+C(2,1))/(4*q(2));
    q(3) = (C(3,2)+C(2,3))/(4*q(2));
    q(4) = (C(3,1)-C(1,3))/(4*q(2));
end
if i == 1 
    q(1) = sqrt(x);
    q(2) = (C(1,2)+C(2,1))/(4*q(1));
    q(3) = (C(1,3)+C(3,1))/(4*q(1));
    q(4) = (C(2,3)-C(3,2))/(4*q(1));
end