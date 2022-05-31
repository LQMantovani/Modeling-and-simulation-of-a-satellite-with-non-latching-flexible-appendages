function [Ml, Ql] = eq_sistem(omega, rot, R, B, Bodies, R_p, U)
%//////////////////////////////////////////////////////////////////%
%System's equations - Generates both generalized mass matrix and 
%generalized forces
%INPUTS:
%omega: Bodies angular velocity in the body reference frame (BRF)
%rot: Attitude quaternions
%R: Bodies inertial position
%B: Elastic coordinates
%Bodies: Data structure with bodies data
%R_p: Bodies velocity in the inertial system (ECI)
%U: Data structure with control forces and torques
%OUTPUTS (non-embedded technique):
%Ml: System's generalized mass matrix
%Ql: System's generalized force vector
%OUTPUTS (embedded technique):
%Ml: States time derivative
%Ql: Empty array
%//////////////////////////////////////////////////////////////////%

muT = Bodies.System.muT;
nflex = Bodies.System.nflex;

% angle2dcm(r1,r2,r3) - rotation order: r1,r2,r3. r1=z, r2=y, r3=x.
C321 = zeros(3,3,length(nflex)); C32 = zeros(3,3,length(nflex)); ArotI_B = zeros(3,3,length(nflex)); 
ArotB_I = zeros(3,3,length(nflex)); g_B = zeros(length(nflex),1); n_B = zeros(length(nflex),1);

%Determines the bodies inertial position, rotation matrices and gravity
for i = 1:length(nflex)
    if Bodies.System.quat == 0 %If Euler angles are used
        phi = rot(1+3*(i-1)); theta = rot(2+3*(i-1)); psi = rot(3+3*(i-1));
        C321(:,:,i) = angle2dcm(psi,theta,phi);
    else %If quaternions are used
        roti = rot(1+4*(i-1):4+4*(i-1));
        Sq = [0 -roti(3) roti(2); roti(3) 0 -roti(1); -roti(2) roti(1) 0];
        C321(:,:,i) = (roti(4)^2-roti(1:3)'*roti(1:3))*eye(3)+2*roti(1:3)*roti(1:3)'-2*roti(4)*Sq;
    end
    %Position variables
    x = R(1,i); y = R(2,i); z = R(3,i); rho = (x^2+y^2+z^2)^.5;
    %Gravity and orbital frequency
    raio = Bodies.System.Radius;
    g_B(i) = muT/raio^2; n_B(i) = -1*sqrt(muT/rho^3);
    %Spherical coordinate angles regarding the LVLH position in ECI system 
    delta = pi/2-acos(z/rho); lambda = atan2(y,x);
    ang = -delta-pi/2;
    C32(:,:,i) = [(-1).*sin(lambda),cos(lambda),0
        (-1).*cos(ang).*cos(lambda),(-1).*cos(ang).*sin(lambda),sin(ang)
        cos(lambda).*sin(ang),sin(ang).*sin(lambda),cos(ang)]; %De forma matricial para deixar a conta mais rapida
    ArotI_B(:,:,i) = C321(:,:,i)*C32(:,:,i);
    ArotB_I(:,:,i) = ArotI_B(:,:,i)';
end

%Mechanisms torques
U_mola_vec = zeros(3,length(nflex),2);
%Verifies if the constraint lock is active or not
if ~Bodies.System.Lock 
    for i = 2:length(nflex)
        if Bodies.B(i).Constraint_type == 1
            %Spring torque acting in the mechanism
            v_i = zeros(3,1); v_i(Bodies.B(i).Body_rotation) = 1;
            v_i = ArotI_B(:,:,1)*ArotI_B(:,:,i)'*v_i;
            %If rotation axis are in oposite directions, the signal must be changed to obtain the correct rotation direction
            signal = v_i(Bodies.B(i).Main_rotation);                %Obtain 1 or -1 to identify if the rotation is positive or negative
            
            Theta_0 = Bodies.B(i).Angular_pos;
            Theta_ref = zeros(3,1); Theta_ref(Bodies.B(i).Body_rotation) = Theta_0(Bodies.B(i).Body_rotation);
            Theta_0 = Theta_0 - Theta_ref;                          %Position the refence vector in a specific position to serve as a reference
            C321_i = angle2dcm(Theta_0(3),Theta_0(2),Theta_0(1))';
            v_id = zeros(3,1); v_id(Bodies.B(i).Axis_Length) = 1;
            v_1 = C321_i*v_id;                                      %Obtains the reference vector in the Body 1 Reference Frame
            v_ii = ArotI_B(:,:,1)*ArotB_I(:,:,i)*v_id;              %Obtains a second vector (aligned with the boom length) in the Body 1 reference frame
            d_Ang_s = acos(dot(v_1, v_ii)/(norm(v_1)*norm(v_ii)));  %Obtains the angle between the vectors

            domega = -(omega(Bodies.B(i).Main_rotation,1)-omega(Bodies.B(i).Body_rotation,i)*signal); %Relative angular velocity in BRF^i
            %Loading mechanism data from data structure
            Theta = Bodies.B(i).Limit_angle;                      %Maximum angle obtained during the experiments
            Theta_r = Bodies.B(i).Angle_rest;                     %Stedy state position obtianed during the expertiments
            prmt = Bodies.B(i).prmt;                              %Model parameters for each mechanism
            d_ang = Theta_r-d_Ang_s;
            
            %First spring behavior (real spring)
            f1 = (prmt(1)/Theta_r)*(Theta-d_Ang_s)*Bodies.B(i).Rot_Dir;
            
            %Second spring behavior (virtual spring)
            f2 = -Bodies.B(i).Rot_Dir*abs((prmt(3)*(Theta_r-d_Ang_s)*sin(prmt(6)*d_Ang_s)))*heaviside(-d_ang);
            
            %Damping behavior
            dmp = -domega*prmt(2)-domega*prmt(5)*(prmt(4)-abs(d_ang))*heaviside(prmt(4)-abs(d_ang));
            
            %TOtatl torque acting is the sum of all three equations above
            U_mola_i = f1 + f2 + dmp;
            
            U_mola_vec(Bodies.B(i).Body_rotation,i,1) = U_mola_i*signal;
            U_mola_vec(Bodies.B(i).Main_rotation,i,2) = -U_mola_i;
            U_mola_vec(:,1,1) = U_mola_vec(:,1,1) + U_mola_vec(:,i,2);
        end
    end
end


U_mola_vec(:,1,1) = U_mola_vec(:,1,1) + U.torque;
F_HUB = U.force; %Appies a force at the HUB

M = []; Q = []; Cq = []; Qc = [];

Q1 = quaternion(ArotB_I(:,:,1)');

if ~Bodies.System.Embedded
    for i = 1:length(nflex)
        %Matrices for each body
        [M_i, Q_i, Cq_i, Qc_i] = Generalized_matrices(i, Bodies, omega, B, C321, C32, ArotB_I, n_B, g_B, U_mola_vec, Q1, R, R_p, F_HUB);
        
        M_temp = zeros(size(M,1),size(M_i,2));
        M = [M M_temp;M_temp' M_i];     %Generalized mass matrix
        Q = [Q;Q_i];                    %Generalized external forces
        Cq = [Cq;Cq_i];                 %Jacobian of the constraint equations
        Qc = [Qc;Qc_i];                 %Right side term associated with the constraints
        
    end
    
    %Output matrix and vector
    Ml = [M Cq';Cq zeros(size(Cq,1))];
    Ql = [Q;Qc];
else
    %Pre-allocating matrices used in Embedded technique
    Cq1 = []; Cq2 = [];
    for i = 1:length(nflex)
        
        [M_i, Qc_i, Q_i, Cq1_i, Cq2_i] = Generalized_matrices_Embedded(i, Bodies, omega, B, C321, C32, ArotB_I, n_B, g_B, U_mola_vec, Q1, R, R_p);
        
        M = blkdiag(M, M_i);
        Q = [Q; Q_i];
        Qc = [Qc; Qc_i];
        Cq1 = [Cq1; Cq1_i];
        Cq2 = blkdiag(Cq2, Cq2_i);
    end
    
    %Output matrix and vector for Embedded technique
    %For Embedded technique, see: Shabana, A. A. Dynamics of multibody systems. ISBN: 9781107337213
    Bdi = [eye(6); -Cq2\Cq1];
    Qcb = [zeros(6,1); Cq2\Qc];
    Mii = Bdi'*M*Bdi;
    Qi  = Bdi'*(Q) - Bdi'*M*Qcb;
    
    q_i = Mii\Qi;                   %Independent coordinates mass matrix
    q_d = -Cq2\Cq1*q_i + Cq2\Qc;    %Dependent coordinates mass matrix
    q = [q_i;q_d];                  %Full system mass matrix
    Ml = q;                         %Returns the stated vector
    Ql = [];
end