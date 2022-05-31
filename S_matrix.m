function [S,A_r,lambda,sigma] = S_matrix(eixo, eixo2, Dims, m, m_tip, I_tt_tip, flag)
% [S,A_r,lambda,sigma] = S_matrix(Axis, Axis2, Dims, m, m_tip, I_tt_tip)
% Function to produce the S matriz based on the modal shapes of the body
%
% INPUT:
% Axis: Axis at which the modal shape is written (between 1 and 3)
% Axis2: Axis that the modal shape infleunces. Must be in format [n1_x, n1_y , n2_z]
% m: mass of the body (without tip mass)
% m_tip: mass of the tip
% I_tt_tip: Inertia of the tip
% OUTPUT:
% S: Modal shape matrix
% A_r: Normalization constants associated with the modal shape 
% lambda: Constants associated with the modal shape
% sigma: Constants associated with the modal shape
%---------------------------------------------%
% Function created in 01/11/2020
% Created by Lorenzzo Mantovani
%---------------------------------------------%

%Verify ir the system has erros;
if eixo<0 || eixo>3
   error('Axis at which the modal shape is written must be between 1 and 3.'); 
end
if length(eixo2) ~= 3
   error('Axis2 vector must have dimension 3.'); 
end
if eixo2(1)< 0 || eixo2(2) < 0 || sum(eixo2()) == 0
   error('Axis2 inputs must be positive integers and its sum must be greater than zero');
end
if eixo2(eixo) ~= 0
   error('Influence of the modal shape in its own axis must be zero (Axis2(Axis) = 0).'); 
end
% if Dims(1)<0||Dims(2)<0||Dims(3)<0||sum(Dims)==0
%     error('Body dimensions must be positive and greater than zero. ');
% end
if m<=0
   error('Mass of the body must be greater than zero.'); 
end
if m_tip<0
   error('m_tip>=0.'); 
end

%Defining the vector of constants for each modal shape
n = max(eixo2);
L_x = Dims(1);
L_y = Dims(2);
L_z = Dims(3);

if flag
    Length = L_x; 
else
    Length = Dims(eixo); 
end


%Position variables
syms x y z zeta

if m_tip ~= 0
    %If there is mass at the tip, uses this solution.
    [~, lambda, sigma] = constants(m/Length, m_tip, I_tt_tip, n, Length);
    phi_r = -(cos(lambda.*zeta)-cosh(lambda.*zeta)+sigma.*(sin(lambda.*zeta)-sinh(lambda.*zeta)));
else
    %If there are not a mass at the tip uses the standard Modal Shapes with
    %the already defined coefficients lambda and sigma.
    %Values obtained from: Blevins, R. D. Formulas for natural frequency and mode shape. DOI: 0442207107
    Lambda_str = [1.875104068710 4.69409113297 7.85475743823 10.9955407348 14.1371683910 17.2787595320];
    Sigma_str  = [.73409551375 1.01846731875 0.99922449651 1.00003355325 0.99999855010 1.00000006265];
    lambda = Lambda_str(1:n)'; sigma = Sigma_str(1:n)';
    phi_r = cosh(lambda.*zeta) - cos(lambda.*zeta) - sigma.*(sinh(lambda.*zeta) - sin(lambda.*zeta));
end

%Obtain the constants in order to make the modal shapes orthogonal
A_r = constants_AR(phi_r, m, m_tip);
phi_r = A_r.*phi_r;
phi_r = transpose(phi_r);    

%S matrix with the contribution of modal shapes in x axis to the y axis
%The partial derivatives terms are multiplied by zero so the Stiffness
%matrix calculated in 'matrices.m' is correct. If the other Stiffness
%matrix is prefered (as presented in Shabana), these multiplications by
%zero must be removed.
if eixo == 1
    phi_r = subs(phi_r,zeta,x/Length);
    if eixo2(2) ~= 0
        S_2 = [-y*diff(phi_r(1:eixo2(2)),x)*0;
            phi_r(1:eixo2(2));
            zeros(1,eixo2(2))];
    else
        S_2 = [];
    end
    if eixo2(3) ~= 0
        S_3 = [-z*diff(phi_r(1:eixo2(3)),x)*0;
            zeros(1,eixo2(3));
            phi_r(1:eixo2(3))];
    else
        S_3 = [];
    end
end

if eixo == 2
    phi_r = subs(phi_r,zeta,y/Length);
    if eixo2(1) ~= 0
        S_2 = [phi_r(1:eixo2(1));
            -x*diff(phi_r(1:eixo2(1)),y)*0;
            zeros(1,eixo2(1))];
    else
        S_2 = [];
    end
    if eixo2(3) ~= 0
        S_3 = [zeros(1,eixo2(3));
            -z*diff(phi_r(1:eixo2(3)),y)*0;
            phi_r(1:eixo2(3))];
    else
        S_3 = [];
    end
end


if eixo == 3
    phi_r = subs(phi_r,zeta,z/Length);
    if eixo2(1) ~= 0
        S_2 = [phi_r(1:eixo2(1));
            zeros(1,eixo2(1));
            -x*diff(phi_r(1:eixo2(1)),z)*0];
    else
        S_2 = [];
    end
    if eixo2(2) ~= 0
        S_3 = [zeros(1,eixo2(3));
            phi_r(1:eixo2(3));
            -y*diff(phi_r(1:eixo2(3)),z)*0];
    else
        S_3 = [];
    end
end

S = [S_2, S_3];