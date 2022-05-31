function [A_r, lambda, sigma] = constants(m, Mt, It, n, L)
%[A_r, lambda, sigma] = constants(m, m_tip, I_tt_tip, n, L)
% Function to find lambdas e sigmas for the beam model used. This function considers the cantilevered-free beam with tip mass
% 
% This function is based on the modal shape discussed by: 
% Erturk, A. Inman, D. J. Piezoelectric Energy Harvesting, 2011. ISBN: 9781119991151
% 
% INPUT:
% m: mass of the body per unit length (withtout the tip mass)
% m_tip: mass of the tip of the body
% I_tt_tip: inertia of the tip
% n: number of flexible modes
% L: length of the body
% OUTPUT:
% A_r: Normalization constants associated with the modal shapes
% lambda: Constants associated with the modal shapes
% sigma: Constants associates with the modal shapes
%---------------------------------------------%
%Function created in 01/11/2020
%Created by Lorenzzo Mantovani
%---------------------------------------------%

%Lambda is determined from the detm = 0. Since mathematica found no
%analytical solution and it depends on the two bodies proprieties, it
%is solved in the beggining of each simulation and stored in Bodies
%structure
detm = @(lm) (-1).*cos(lm).^2+(-1).*It.*L.^(-4).*lm.^4.*m.^(-2).*Mt.*cos( ...
    lm).^2+(-2).*cos(lm).*cosh(lm)+2.*It.*L.^(-4).*lm.^4.*m.^( ...
    -2).*Mt.*cos(lm).*cosh(lm)+(-1).*cosh(lm).^2+(-1).*It.*L.^( ...
    -4).*lm.^4.*m.^(-2).*Mt.*cosh(lm).^2+2.*It.*L.^(-3).*lm.^3.* ...
    m.^(-1).*cosh(lm).*sin(lm)+2.*L.^(-1).*lm.*m.^(-1).*Mt.* ...
    cosh(lm).*sin(lm)+(-1).*sin(lm).^2+(-1).*It.*L.^(-4).* ...
    lm.^4.*m.^(-2).*Mt.*sin(lm).^2+2.*It.*L.^(-3).*lm.^3.*m.^( ...
    -1).*cos(lm).*sinh(lm)+(-2).*L.^(-1).*lm.*m.^(-1).*Mt.*cos( ...
    lm).*sinh(lm)+sinh(lm).^2+It.*L.^(-4).*lm.^4.*m.^(-2).*Mt.* ...
    sinh(lm).^2;

%Initializing the vectors
lambda = zeros(1,n); sigma = lambda; AC = zeros(size(lambda,2),2);

i = 1;
c0 = .1;
while i <= size(lambda,2)
    lambda(i) = fsolve(detm,c0,optimset('Disp','off')); %Solving lambda
    if i == 1
        if lambda(i) ~= 0
            i=i+1; %If lambda is solved
        end
    else
        if lambda(i) > lambda(i-1)*1.1
            i = i+1;
        end
    end
    c0 = c0+.1; %Increase the intial guess to the next solution be found
end
for i = 1:size(lambda,2)
    
    %Solving for sigma
    sigma(i) = (sin(lambda(i))-sinh(lambda(i))+lambda(i)*(Mt/(m*L))*(cos(lambda(i))-cosh(lambda(i))))...
        /(cos(lambda(i))+cosh(lambda(i))-lambda(i)*(Mt/(m*L))*(sin(lambda(i))-sinh(lambda(i))));
    
    %Solving for the AC constants
    lm = lambda(i);
    matrix = [cos(lm)+cosh(lm)+(-1).*It.*L.^(-3).*lm.^3.*m.^(-1).*(sin( ...
        lm)+sinh(lm)),It.*L.^(-3).*lm.^3.*m.^(-1).*(cos(lm)+(-1).* ...
        cosh(lm))+sin(lm)+sinh(lm);L.^(-1).*lm.*m.^(-1).*Mt.*(cos( ...
        lm)+(-1).*cosh(lm))+sin(lm)+(-1).*sinh(lm),(-1).*cos(lm)+( ...
        -1).*cosh(lm)+L.^(-1).*lm.*m.^(-1).*Mt.*(sin(lm)+(-1).*sinh( ...
        lm))];
    
    Modal_shapes_AC = @(A) [A(1)*matrix(1,1)+A(2)*matrix(1,2); A(1)*matrix(2,1)+A(2)*matrix(2,2)];
    AC(i,:) = fsolve(Modal_shapes_AC,[1 1],optimset('Disp','off'));
    A_r = AC(:,1);
end
lambda = lambda';
sigma = sigma';