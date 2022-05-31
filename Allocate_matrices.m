%---------------------------------------%
% For each body with flexible characteristics, generates the mass matrices
% and saves it in the Data structure Bodies in use.
% First verifies if the existing Data structure Bodies already has the
% matrices to avoid unecessary calculations. At the end saves an updates
% Bodies Data Structure.
%
%---------------------------------------------%

% Verify if the existing Bodies Data Structure already has the necessary mass matrices
flag_file = false;
directory = dir;
for i=1:size(directory,1)
    if directory(i).name == "Bodies_data.mat"
        flag_file = true;
        load('Bodies_data.mat');
    end 
end

% In case there is a saved file, compares the fields data to see is
% fomething changed. If nothing changed, the previous file can be used.
Change = 1:length(nflex);
flag_change(1:length(nflex)) = false;
if flag_file 
    for i = 2:min([length(Bodies_stored.B) length(nflex)])
        if nflex(i)
            if ~all(Bodies_stored.B(i).Flexible_modes == Bodies.B(i).Flexible_modes)
                flag_change(i) = true;
            end
            if Bodies_stored.B(i).mass ~= Bodies.B(i).mass
                flag_change(i) = true;
            end
            if Bodies_stored.B(i).Tip_mass ~= Bodies.B(i).Tip_mass
                flag_change(i) = true;
            end
            if Bodies_stored.B(i).Tip_Inertia ~= Bodies.B(i).Tip_Inertia
                flag_change(i) = true;
            end
            Dim = [Bodies.B(i).Length Bodies.B(i).Internal_Radius Bodies.B(i).External_Radius];
            Dim_stored = [Bodies_stored.B(i).Length Bodies_stored.B(i).Internal_Radius Bodies_stored.B(i).External_Radius];
            if Dim(1) ~= Dim_stored(1) || Dim(2) ~= Dim_stored(2) || Dim(3) ~= Dim_stored(3)
                flag_change(i) = true;
            end
            if Bodies_stored.B(i).Axis_Length ~= Bodies.B(i).Axis_Length
                flag_change(i) = true;
            end
            if Bodies_stored.B(i).Flexible_modes ~= Bodies.B(i).Flexible_modes
                flag_change(i) = true;
            end
            if Bodies_stored.B(i).Young_module ~= Bodies.B(i).Young_module
                flag_change(i) = true;
            end
            if Bodies_stored.B(i).Poisson ~= Bodies.B(i).Poisson
                flag_change(i) = true;
            end
            if ~flag_change(i)
                %If there are no changes, copies the Func from the stored Bodies to the current Bodies.
                Bodies.B(i).Func = Bodies_stored.B(i).Func;
                if length(Bodies.B(i).Structural_Damping) == 1
                    %Modal damping
                    Bodies.B(i).Func.D_ff = double(Bodies.B(i).Structural_Damping*eye(size(Bodies.B(i).Func.K_ff, 1)));
                else
                    %Proportional damping
                    Bodies.B(i).Func.D_ff = double(Bodies.B(i).Structural_Damping(1)*Bodies.B(i).Func.K_ff + Bodies.B(i).Structural_Damping(2)*Bodies.B(i).Func.M_FF);
                end
            end
        end
    end
    for i = (length(Bodies_stored.B)+1):length(nflex)
        flag_change(i) = true;
    end
else
   flag_change(1) = false; %standard so the first body is not evaluated
   flag_change(2:length(nflex)) = true; 
end

%If there are changes in any of the Bodies, it generates the Generalized
%Matrices.
for i = Change(flag_change)
    if Bodies.B(i).nflex
        disp('---------------------------------------------------------')
        disp(['Updating Generalized Matrices for Body ',num2str(i),' .']);
        m        = Bodies.B(i).mass;
        m_tip    = Bodies.B(i).Tip_mass;
        I_TT_tip = Bodies.B(i).Tip_Inertia;
        Dim      = [Bodies.B(i).Length Bodies.B(i).Internal_Radius Bodies.B(i).External_Radius];
        Axis     = Bodies.B(i).Axis_Length;
        Axis2    = Bodies.B(i).Flexible_modes;
        Young    = Bodies.B(i).Young_module;
        Poisson  = Bodies.B(i).Poisson;
        [M_RR, I_TT, I_TF, M_FF, S_B, S_TB_SK, Qvr, Qva, Qvf, K_ff, qf, qf_dot, SM, k, G, CG] = matrices(m, m_tip, I_TT_tip, Dim, Axis, Axis2, Young, Poisson, true, Bodies, i);
        
        syms omega_x omega_y omega_z theta1 theta2 theta3 theta4 x y z 
        vec = sym('vec', [3,1]);
        
        Bodies.B(i).Func.M_RR = double(M_RR);
        Bodies.B(i).Func.I_TT = matlabFunction(I_TT,'Vars',{qf});
        Bodies.B(i).Func.I_TF = matlabFunction(I_TF,'Vars',{qf});
        Bodies.B(i).Func.M_FF = double(M_FF);
        Bodies.B(i).Func.S_B  = double(S_B);
        Bodies.B(i).Func.S_TB_SK = matlabFunction(S_TB_SK,'Vars',{qf});
        Bodies.B(i).Func.Qvr = matlabFunction(Qvr,'Vars',{[omega_x,omega_y,omega_z],qf,qf_dot,[theta1;theta2;theta3;theta4]});
        Bodies.B(i).Func.Qva = matlabFunction(Qva,'Vars',{[omega_x,omega_y,omega_z],qf,qf_dot});
        Bodies.B(i).Func.Qvf = matlabFunction(Qvf,'Vars',{[omega_x,omega_y,omega_z],qf,qf_dot});
        Bodies.B(i).Func.K_ff = double(K_ff);
        Bodies.B(i).Func.D_ff = double(Bodies.B(i).Structural_Damping(1)*K_ff + Bodies.B(i).Structural_Damping(2)*M_FF); %Alpha and Beta parameters must be incliuded in Bodies_Data fields
        Bodies.B(i).Func.S = matlabFunction(SM, 'Vars', {[x, y, z]});
        Bodies.B(i).Func.k = matlabFunction(k, 'Vars', {vec,[theta1; theta2; theta3; theta4]});
        Bodies.B(i).Func.G = matlabFunction(G, 'Vars', {[theta1; theta2; theta3; theta4]});
        Bodies.B(i).Func.Center_of_mass = matlabFunction(CG, 'Vars', {qf});
        
        disp('Update finished');
        disp('---------------------------------------------------------')
    end
end


%At the end, saves Bodies to be used further.
Bodies_stored = Bodies;
save('Bodies_data','Bodies_stored');