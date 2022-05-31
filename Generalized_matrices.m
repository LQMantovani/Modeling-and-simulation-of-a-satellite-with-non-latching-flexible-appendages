function [M, Q, Cq, Qc] = Generalized_matrices(i, Bodies, omega, B, C321, C32, ArotB_I, n_B, g_B, U_mola_vec, Q1, R, Rp, F_HUB)
%//////////////////////////////////////////////////////////////////%
%Function used to generate the mass, inertia and flexible mass
% matrices, generalized forces and constraints
%
%INPUTS:
%i: Body number
%Bodies: Data structure with bodies information
%omega: Body angular velocity in the body fixed reference frame (BRF)
%B: Data structure with elastic coordinates
%C321: Roation matrices from the LVLH frame to the BRF
%C32: Rotation matrices from the ECI to the LVLH frame
%ArotB_I: Rotation matrices from the BRF to the ECI
%n_B: Orbit's frequency
%g_B: Gravity acceleration at the BRF
%U_mola_vec: Matrix with mechanism torques in each body
%Q1: Main body attitue quaternions
%R: Body inertial position
%Rp: Body inertial velocity
%F_HUB: Force actuating in the HUB
%
%OUTPUT:
%M: Generalized mass matrix
%Q: Generalized forces vector
%Cq: Jacobian of the contraint equations
%Qc: "Force" vector associated with the constraints
%//////////////////////////////////////////////////////////////////%

nflex = Bodies.System.nflex;
Qi = quaternion(ArotB_I(:,:,i)');

if Bodies.B(i).nflex == 0
    
    %For rigid bodies with reference frames at their centers of mass
    %Matrix of translational inertia
    mrr = [1 0 0;0 1 0;0 0 1]*Bodies.B(i).mass;
    %Matrix of rotational inertia
    Itt = Bodies.B(i).Inertia;
    %General Mass matrix
    M = [mrr zeros(3);zeros(3) Itt];
    
    %FORCES
    Qer = C32(:,:,i)'*[0; 0; Bodies.B(i).mass*g_B(i)];
    %TORQUES
    Qva = -cross(omega(:,i),Itt*omega(:,i));
    %Gravitaional torque
    %See: Wie, Bong. Space Vehicle Dynamics and Control, second edition.
    %DOI: 10.2514/4.860119
    Ra = [0;0;-1];
    Rb = C321(:,:,i)*Ra;
    Qg = 3*n_B(i)^2*cross(Rb,Itt*Rb);
    
    Qea=U_mola_vec(:,i,1);
    
    if i ~= 1, F_HUB = F_HUB*0; end %Applies a force only at the HUB
    
    Q=[Qer + C32(:,:,i)'*C321(:,:,i)'*F_HUB
        Qva + Qg + Qea];
    
else
    %For flexible bodies with generic body fixed frame origin
    qf = B(i).Qf;
    qfp = B(i).Qfp;
    mrr = Bodies.B(i).Func.M_RR;
    Sttx = Bodies.B(i).Func.S_TB_SK(qf);
    St =  Bodies.B(i).Func.S_B;
    Itt = Bodies.B(i).Func.I_TT(qf);
    Itf = Bodies.B(i).Func.I_TF(qf);
    mff = Bodies.B(i).Func.M_FF;
    Qvr = Bodies.B(i).Func.Qvr(omega(:,i)', qf, qfp, Qi');
    Qva = Bodies.B(i).Func.Qva(omega(:,i)', qf, qfp);
    Qvf = Bodies.B(i).Func.Qvf(omega(:,i)', qf, qfp);
    Kff = Bodies.B(i).Func.K_ff;
    Dff = Bodies.B(i).Func.D_ff;
    
    %Generalized mass matrix
    M=[mrr (ArotB_I(:,:,i)*Sttx') ArotB_I(:,:,i)*St
        (ArotB_I(:,:,i)*Sttx')'   Itt  Itf
        (ArotB_I(:,:,i)*St)'  Itf' mff];
    
    %FORCES
    %Gravity force in the body system
    Fg = C321(:,:,i)*[0; 0; (Bodies.B(i).mass + Bodies.B(i).Tip_mass)*g_B(i)];
    
    %Generalized force associated with gravity
    % In accordance with Shabana - Dynamics of Multibody Systems (2013). Example 5.5
    Qrg = ArotB_I(:,:,i)*Fg;
    
    %Genralized forces in the inertial frame
    Qtg = -(Fg'*Sttx)'/(Bodies.B(i).mass + Bodies.B(i).Tip_mass);
    %Influence of gravity on flexible modes
    Qfg = (Fg'*St)'/(Bodies.B(i).mass + Bodies.B(i).Tip_mass);
    
    %TORQUES
    %Gravitaional torque
    Ra = [0;0;-1];
    Rb = C321(:,:,i)*Ra;
    Qg = 3*n_B(i)^2*cross(Rb,Itt*Rb);
    
    %Generalized external forces
    Qer = zeros(3,1);
    Qea = U_mola_vec(:,i,1);
    
    Q=[Qvr + Qrg + Qer
        Qva + Qtg + Qg + Qea
        Qvf + Qfg - Kff*qf - Dff*qfp];
    
end


%Generate the constraints - except for the first body (it is not required)
if i ~= 1
    if Bodies.B(i).Constraint_type == 1
        
        %Pre-allocate the matrix for performance improvment
        Cq = zeros(5,length(nflex)*6+sum(nflex));
        
        %Obtain the quaternions outside the constraint functions to avoid
        %unecessary calls and improve code performance
        
        %Obtains the parameters to constraint three translational degrees
        %of freedom
        %See: Shabana, A. A. Computational Dynamics. ISBN: 9780470686850. DOI: 10.1002/9780470686850
        [Gb_1, k_1, Cqf_1] = constraints_trl_quat(omega(:,1), Q1, Bodies.B(i).Body_position);
        [Gb_i, k_i, Cqf_i] = constraints_trl_quat(omega(:,i), Qi, Bodies.B(i).Main_position);
        
        %Body states position in the X vector
        tr_i = Bodies.B(i).pos_v(1);
        rot_i = Bodies.B(i).pos_omega(1);
        
        %Translation constraints
        Cq(1:3,1:3) = -[1 0 0;0 1 0;0 0 1];   %Always in this position, since it is the constraint associated with the first body
        Cq(1:3,4:6) = -k_1*Gb_1;              %Always in this position, since it is the constraint associated with the first body
        Cq(1:3,tr_i:tr_i+2) = [1 0 0;0 1 0;0 0 1];
        Cq(1:3,rot_i:rot_i+2) = k_i*Gb_i;
        
        %Translation constraint forces - The same for locked and free mechanisms
        %Adds the contraints control to mitigate the contraints drift over time
        %See: Shabana, A. A. Computational Dynamics. ISBN: 9780470686850. DOI: 10.1002/9780470686850
        %and: Bauchau, O.; Laulusa, A. Review of Contemporary Approaches. DOI: 10.1115/1.2803258
        c = R(:,1) + ArotB_I(:,:,1)*Bodies.B(i).Body_position - (R(:,i) + ArotB_I(:,:,i)*Bodies.B(i).Main_position);
        cp = Rp(:,1) + k_1*Gb_1*omega(:,1) - (Rp(:,i) + k_i*Gb_i*omega(:,i));
        K = Bodies.System.K_p;
        Qc(1:3,1) = -Cqf_i-(-Cqf_1) + 2*K*cp + K^2*c;
        
        %Rotation constraint forces
        Kr = Bodies.System.K_r; %Constraint drift control gain
        if ~Bodies.System.Lock
            %In case the body is not locked
            %See: Shabana, A. A. Computational Dynamics. ISBN: 9780470686850. DOI: 10.1002/9780470686850
            [Cq_1, Cq_i, Cqf] = constraints_rot_quat(omega(:,1), omega(:,i), Q1, Qi, Bodies.B(i).Body_rotation, Bodies.B(i).Main_rotation);
            
            crp = Cq_1*omega(:,1) + Cq_i*omega(:,i);
            
            Cq(4:5,4:6) = Cq_1;               %Always in this position, since it is the constraint associated with the first body
            Cq(4:5,rot_i:rot_i+2) = Cq_i;
            
            %Constraints contol to mitigate drift
            va = [0;0;0]; va(Bodies.B(i).Body_rotation) = 1; 
            
            switch Bodies.B(i).Main_rotation
                case 1
                    vb = [0; 1; 0]; vc = [0; 0; 1];
                case 2
                    vb = [1; 0; 0]; vc = [0; 0; 1];
                case 3
                    vb = [1; 0; 0]; vc = [0; 1; 0];
            end
            
            %Constraints contol to mitigate drift
            cr = [dot(ArotB_I(:,:,1)'*ArotB_I(:,:,i)*va, vb); dot(ArotB_I(:,:,1)'*ArotB_I(:,:,i)*va, vc)];
            Qc(4:5,1) = Cqf - 2*Kr*crp - Kr^2*cr;
            
        else            
            %For a locked body
            %Defines three orthogonal vectores
            a_j_x = [1;0;0]; a_j_y = [0;1;0]; a_j_z = [0;0;1];
            %Arrange the three orthogonal vectors in a matrix
            a_j = [a_j_x, a_j_y, a_j_z];
            %Define parallel vector in the attached body
            a_i_p = Bodies.B(1).ArotB_I'*Bodies.B(i).ArotB_I*a_j;
            %Rearange the vector in matrix so each collum in a_i is
            %orthogonal to each collum in a_j
            a_i = [a_i_p(:,2), a_i_p(:,3), a_i_p(:,1)];
            %Calculate the matrix terms
            [Cq_i, Cq_j, Cqf] = constraint_fixed(a_i, a_j, Q1', Qi', omega(:,1), omega(:,i));
            Cq(4:6,4:6) = Cq_i;
            Cq(4:6,rot_i:rot_i+2) = Cq_j;
            %Calculate the equations
            cr = [dot(ArotB_I(:,:,1)*a_i(:,1),ArotB_I(:,:,i)*a_j(:,1))
                dot(ArotB_I(:,:,1)*a_i(:,2),ArotB_I(:,:,i)*a_j(:,2))
                dot(ArotB_I(:,:,1)*a_i(:,3),ArotB_I(:,:,i)*a_j(:,3))];
            %Calculate the equation derivative
            crp = Cq_i*omega(:,1) + Cq_j*omega(:,i);
            Qc(4:6) = Cqf - 2*Kr*crp - Kr^2*cr;
        end
    end
else
    Cq = []; Qc = [];
end