function I = int_v(F, dx, dy, dz, Lx, Ly, Lz, flag, Axis)
% I = int_v(F, dx, dy, dz, Lx, Ly, Lz)
% Function to perform triple integration
% INPUT:
% F: Function to integrate
% dx: first variable to integrate over
% dy: second variable to integrate over
% dz: third variable to integrate over
% Lx: First integration limit [Lx_inf Lx_sup]
% Ly: First integration limit [Ly_inf Ly_sup]
% Lz: First integration limit [Lz_inf Lz_sup]
% flag = true if cilindrical coordinates should be used
% OUTPUT:
% I: Integrated value
%---------------------------------------------%
%Function created in 01/11/2020
%Created by Lorenzzo Mantovani
%---------------------------------------------%

if ~flag
    %In case of cartesian integration
    P1 = int(F,dx,Lx);
    P2 = int(P1,dy,Ly);
    I = real(vpa(simplify(int(P2,dz,Lz))));
    
else
    %In case of cilindrical coordinates integration
    %Lists the variables available to pick the correct one to perform the integration
    Vars = [dx dy dz];
    %Defines the length
    L = Lx;
    %Defines the internal and external radius
    r_i = Ly; r_e = Lz;
    if isnumeric(r_i) && r_i < 0
        error('Internal radius must be equal or greater than zero.')
    end
    
    syms raio phi
    
    switch Axis
        case 1
            %Caso o comprimento seja ao longo do eixo X:
            F_t1 = subs(F, dy, raio*cos(phi));
            F_sph = raio*subs(F_t1, dz, raio*sin(phi));
        case 2
            F_t1 = subs(F, dx, raio*cos(phi));
            F_sph = raio*subs(F_t1, dz, raio*sin(phi));
        case 3
            F_t1 = subs(F, dx, raio*cos(phi));
            F_sph = raio*subs(F_t1, dy, raio*sin(phi));
        otherwise
            error('Invalid Axis value. Axis value must be in range [1 3].');
    end
    
    P1 = int(vpa(F_sph), Vars(Axis),[0 L]);
    P2 = int(P1, raio, [r_i r_e]);
    I  = real(vpa(simplify(int(P2, phi, [0 2*pi]))));
end