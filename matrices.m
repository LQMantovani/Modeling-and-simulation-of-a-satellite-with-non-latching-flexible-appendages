function [M_RR, I_TT, I_TF, M_FF, S_B, S_TB_SK, Qvr, Qva, Qvf, K_ff, qf, qf_dot, SM, k, G, CG] = matrices(m, m_tip, I_TT_tip, Dim, Axis, Axis2, Young, Poisson, flag, Bodies, Body_number)
%[I_tt, I_tf, m_ff, S_b, S_tb_sk,Qvr,Qva,Qvf,K_ff,qf,qf_dot] = matrices(m, m_tip, I_tt_tip, Dim, Axis, Axis2)
%Code to obtain the mass matrices of a generic body with flexible modes
%Flexible modes are of a cantilevered beam with tip mass and inertia
% INPUT:
% m: mass of the body (without tip mass)
% m_tip: mass of the tip
% I_tt_tip: Inertia of the tip (scalar, is assumed the same for all axis)
% Dim: Dimensions of the body [X, Y, Z]
% Axis: Axis at which the modal shape is written (between 1 and 3)
% Axis2: Axis that the modal shape infleunces. Must be in format [n_x, n_y , n_z]
% Young: Young's Modulus
% Poisson: Poisson
% flag: Set to true if the structure is cilindric and cilindric coordinates shoud be used. Note that 
%    in this case Dim(1) is the length, Dim(2) is the internal radius and Dim(3) the external Radius
% Bodies: Data structure with bodies data
% Body_number: Number of the current body being evaluated
% OUTPUT:
% M_RR: Translational mass matrix
% I_TT: Inertia matrix
% I_TF: Coupling between rotation and flexible mass matrix
% M_FF: Modal mass matrix
% S_B: Couplying between transaltion and flexible modes
% S_TB_SK: Couplying between translation and rotation matrix
% Qvr: Translational generalized forces
% Qva: Rotational generalized forces
% Qvf: Flexible modes generalized forces
% K_ff: Sitffness matrix
% qf: Generalized coordinates vector. Return the qf vector as qf = [qf_x; qf_y; qf_z] as standard.
% qf_dot: Generalized accelerations vector
% SM: Modal Shape matrix
% k: Jacobian matrix
% G: Matrix relating angular velocity and attitude parameters time
%   derivative
% CG: Center of mass position vector
%---------------------------------------------%
%Function created in 01/11/2020
%Created by Lorenzzo Mantovani
%---------------------------------------------%


%Verify ir the system has erros;
if Axis<1 || Axis>3
   error('Axis at which the modal shape is written must be between 1 and 3.'); 
end
if length(Axis2) ~= 3
   error('Axis2 vector must have dimension 3.'); 
end
if Axis2(1)< 0 || Axis2(2) < 0 || sum(Axis2()) == 0
   error('Axis2 inputs must be positive integers and its sum must be greater than zero');
end
if Axis2(Axis) ~= 0
   error('Influence of the modal shape in its own axis must be zero (Axis2(Axis) = 0).'); 
end
if isnumeric(Dim) && Dim(1)<0||Dim(2)<0||Dim(3)<0||sum(Dim)==0
    error('Body dimensions must be positive and greater than zero. ');
end
if m <= 0
   error('Mass of the body must be greater than zero.'); 
end
if m_tip < 0
   error('m_tip >= 0.'); 
end


%Matrices implemented here are based on: 
% Shabana, A. A. Dynamics of Multibody Systems. 2013. ISBN: 9781107337213

%Number of flexible modes desired:
n = sum(Axis2);
L_x = Dim(1);
L_y = Dim(2);
L_z = Dim(3);

%Position coordinates in the inertial frame
syms X Y Z

R = [X; Y; Z];

%Underformed position of the point P in Body i in the BRF
syms x y z

u_0 = [x; y; z];

%Number of generalized elastic coordinates
qf = sym('qf',[n,1]);

%S matrix of modal shapes
S = sym('S',[3,n]);

%Quaternion parametrization
theta = sym('theta',[4,1],'real');

%Rotation matriz from the inertial frame to the BRF
Sq=[0 -theta(3) theta(2); theta(3) 0 -theta(1); -theta(2) theta(1) 0];
AT=(theta(4)^2-theta(1:3)'*theta(1:3))*eye(3)+2*theta(1:3)*theta(1:3)'-2*theta(4)*Sq;
        
%Rotation matriz from the BRF to the inertial frame
A = transpose(AT);

%Position of point P in Body i in the inertial frame
u = u_0 + S*qf;

%%
%Defining the generalized coordinates of the Body i
q = [R; theta; qf];

%Time derivatives of the coordinates
syms X_dot Y_dot Z_dot x_dot y_dot z_dot

R_dot = [X_dot; Y_dot; Z_dot];
theta_dot = sym('theta_dot',[4,1]);
qf_dot = sym('qf_dot',[n,1]);
% u_0_dot = [x_dot; y_dot; z_dot];

q_dot = [R_dot; theta_dot; qf_dot];

%Matriz used to obtain genrealized forces
vec = sym('vec', [3,1]);
k = jacobian(A*vec, theta);
G = 1/2*[theta(4) -theta(3) theta(2)
    theta(3) theta(4) -theta(1)
    -theta(2) theta(1) theta(4)
    -theta(1) -theta(2) -theta(3)];
%%
%Section to obtain the generalized mass matrix

%Spliting the generalized mass matriz in smaller parts
%Only considers the mass influence on the modal shape if the mass is at
%the bodie's tip
if Bodies.B(Body_number).Length == Bodies.B(Body_number).Tip_pos
    [SM,~,~,~] = S_matrix(Axis,Axis2,Dim, m, m_tip, I_TT_tip, flag);
else
    
    [SM,~,~,~] = S_matrix(Axis,Axis2,Dim, m, m_tip*0, I_TT_tip*0, flag);
end


%Defines the integration limits 
if ~flag
    lim_x = [0 L_x];
    lim_y = [0 L_y];
    lim_z = [0 L_z];
    %Density
    rho = m/(L_x*L_y*L_z);
else
    %Defines the length
    Lg = Dim(1);
    %Defines the internal and external radius
    r_i = Dim(2); r_e = Dim(3);
    %Density
    rho = m/(Lg*pi*(r_e^2-r_i^2));
%     rho = Bodies.B(Body_number).Density;
    lim_x = Lg;
    lim_y = r_i;
    lim_z = r_e;
end

%Skew symmetric matrix of the u vector
u = subs(u,S,SM);
u_sk = [0 -u(3) u(2); u(3) 0 -u(1); -u(2) u(1) 0];

%Generalized matrices:
% m_rr = eye(3)*m;
m_rr = eye(3)*int_v(rho, x, y, z, lim_x, lim_y, lim_z, flag, Axis);
%Inertia
I_tt  = int_v(rho*transpose(u_sk)*u_sk,x,y,z, lim_x, lim_y, lim_z, flag, Axis);
%Coupling between inertia and flexible mode
I_tf = int_v(rho*u_sk*SM,x,y,z, lim_x, lim_y, lim_z, flag, Axis);
%Coupling between translation and rotatin
S_tb_sk = int_v(rho*u_sk,x,y,z, lim_x, lim_y, lim_z, flag, Axis);
%Coupling between translation and flexible mode
S_b = int_v(rho*SM,x,y,z, lim_x, lim_y, lim_z, flag, Axis);
%Mass matrix associated with the flexible modes
m_ff = int_v(rho*transpose(SM)*SM,x,y,z, lim_x, lim_y, lim_z, flag, Axis);

S_tb = int_v(rho*u,x,y,z, lim_x, lim_y, lim_z, flag, Axis);
%%
%To add up the influence of the tip mass on the mass matrix.
%Considers an mass point
m_rr_tip = eye(3)*m_tip;

%Position vector of the tip mass
if ~flag
    Lg = Dim(Axis);
end

u_0_tip = zeros(3,1);
u_0_tip(Axis) = Bodies.B(Body_number).Tip_pos;

Limits = zeros(1,3); Limits(Axis) = Bodies.B(Body_number).Tip_pos;

u_tip = u_0_tip + subs(SM*qf,[x,y,z],Limits);

%Skew symmetric matrix of the tip's mass position vector
u_sk_tip = [0 -u_tip(3) u_tip(2); u_tip(3) 0 -u_tip(1); -u_tip(2) u_tip(1) 0];

%Inertia matrix associated with the tip mass
I_tt_tip = m_tip*transpose(u_sk_tip)*u_sk_tip;
%Coupling between inertia and flexible mode for the tip mass
I_tf_tip = m_tip*subs(u_sk_tip*SM,[x,y,z],Limits);

S_tb_tip = m_tip*u_tip;
%Coupling between translation and rotatin
S_tb_sk_tip = m_tip*subs(u_sk_tip,[x,y,z],Limits);
%Coupling between translation and flexible mode
S_b_tip = m_tip*subs(SM,[x,y,z],Limits);
%Mass matrix associated with the flexible modes
m_ff_tip = m_tip*subs(transpose(SM)*SM,[x,y,z],Limits);
%%
%The final mass matrix with tip mass:

M_RR = m_rr + m_rr_tip;
I_TT = I_tt + I_tt_tip ;
I_TF = I_tf + I_tf_tip;
M_FF = m_ff + m_ff_tip;
S_B  = S_b + S_b_tip;
S_TB = S_tb + S_tb_tip;
S_TB_SK = S_tb_sk + S_tb_sk_tip;

%%
%To add the influence of a tip inertia on the system

%Rotation of the modal shape
dS = diff(SM,y);
%Create a vector pointing towards the booms tip to obtain the correct
%rotation axis
vec_length = zeros(3,1); vec_length(Axis) = 1;
vec_length_skew = [0 -vec_length(3) vec_length(2); vec_length(3) 0 -vec_length(1); -vec_length(2) vec_length(1) 0];
%Rotation of the tip of the body
dS_tip = vec_length_skew*subs(dS,[x,y,z], Limits);

I_TT_TIP = eye(3)*I_TT_tip;
I_TF_TIP = I_TT_TIP*dS_tip;
M_FF_TIP = transpose(dS_tip)*I_TT_TIP*dS_tip;

I_TT = I_TT + I_TT_TIP;
I_TF = I_TF + I_TF_TIP;
M_FF = M_FF + M_FF_TIP;
%%
%Force matrices

%Angular velocity
syms omega_x omega_y omega_z
omega = [omega_x; omega_y; omega_z];

for ii = 1:3
    for jj = 1:3
        dI_TT(ii,jj) = jacobian(I_TT(ii,jj),q)*q_dot;
    end
end
omega_sk = [0 -omega(3) omega(2); omega(3) 0 -omega(1); -omega(2) omega(1) 0];
Qvr = -A*(omega_sk*omega_sk*S_TB + 2*omega_sk*S_B*qf_dot);
Qva = -cross(omega,I_TT*omega) - dI_TT*omega - cross(omega,I_TF*qf_dot);
arg = rho*transpose(SM)*(omega_sk*omega_sk*u + 2*omega_sk*SM*qf_dot);
arg_tip = transpose(SM)*(omega_sk*omega_sk*u_tip + 2*omega_sk*SM*qf_dot);
Qvf = -int_v(arg,x,y,z, lim_x, lim_y, lim_z, flag, Axis) - m_tip*subs(arg_tip,[x y z],Limits);

%Stifness matrix
%This method introduces the area moment of inertia automatically, but the
%system becomes stiffer. It acompanies changes in S-matrices (requires the
%axial compression to be taking in account)
% Lambda = Young*Poisson/((1+Poisson)*(1-2*Poisson));
% Mu = Young/(2 + 2*Poisson);
% E = blkdiag([Lambda+2*Mu Lambda Lambda; Lambda Lambda+2*Mu Lambda; Lambda Lambda Lambda+2*Mu],diag(2*Mu*ones(3,1)));
% DS = 0.5*[2*diff(SM(1,:),x); 2*diff(SM(2,:),y); 2*diff(SM(3,:),z); diff(SM(1,:),y) + diff(SM(2,:),x);diff(SM(1,:),z) + diff(SM(3,:),x); diff(SM(2,:),z) + diff(SM(3,:),y)];
% 
% K_ff = int_v(transpose(DS)*E*DS,x,y,z, lim_x, lim_y, lim_z, flag, Axis);

%Then, this method will be used similar to: Murilo, A.; et al. Real-time
%implementation of a parameterized Model Predictive Control for Attitude
%Control Systems of rigid-flexible satellite. DOI: 10.1016/j.ymssp.2020.107129
%Area moment of inertia
I_A = pi/4*(r_e^4 - r_i^4);
%Stiffness matrix - This must be written in a generalized format later
K_ff = Young*I_A*int(transpose(diff(SM,'y',2))*diff(SM,'y',2),y, [0 Lg]);
%%
%In order to avoid some imaginary number than can occur due to numeric
%erros and simplify the results to improve further performance since they
%are called several times during the simulation.

M_RR = real(vpa(M_RR));
I_TT = real(vpa(simplify(I_TT)));
I_TF = real(vpa(simplify(I_TF)));
M_FF = real(vpa(simplify(M_FF)));
S_B = real(vpa(simplify(S_B)));
S_TB_SK = real(vpa(simplify(S_TB_SK)));
Qvr = real(vpa(simplify(Qvr)));
Qva = real(vpa(simplify(Qva)));
Qvf = real(vpa(simplify(Qvf)));
K_ff = real(vpa(simplify(K_ff)));

%%
%To determine the center of mass position
CG = S_TB/M_RR(1,1);