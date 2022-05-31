function AR = constants_AR(phi, m, m_tip)
%AR = constants_AR(phi, m, m_tip)
%-----------------------------------------------------------%
% Function to normilize the modal shape
% 
% This function is based on the mass normalization presented by: 
% Murilo, A.; et al. Real-time implementation of a parameterized Model Predictive Control for 
% Attitude Control Systems of rigid-flexible sattelite. DOI: 10.1016/j.ymssp.2020.107129
% 
% INPUT:
% phi: Modal shapes
% m: mass
% m_tip: tip_mass
% OUTPUT:
% AR: Normalization constant associated with the modal shape
%-----------------------------------------------------------%

%Allocates the vector
AR = sym('AR', [length(phi), 1]);

syms zeta

%Allocate the vector with modal shape and symbolic normalization constants
phi_r = AR.*phi;

aux = 1;
for i = 1:length(AR)
   for j = 1:length(AR)
            eq(aux) = m*int(phi_r(i)*phi_r(j), zeta, [0, 1]) + m_tip*subs(phi_r(i)*phi_r(j), zeta, 1);
        if i == j
            eq(aux) = eq(aux)-(m+m_tip);
        end
        aux = aux + 1;
   end
end

eq2 = matlabFunction(vpa(eq), 'Vars', {AR});

%Obtains the normalization constants
AR = fsolve(eq2, ones(length(AR),1),optimset('Display','off','Diagnostics','off','Algorithm','levenberg-marquardt'));