%Code used to obtain the relative angular position, relative angular
%velocity, relative linear position and its erros and system angular
%momentum.

%Initializing variables
P2_1 = zeros(length(t),3,length(nflex)); P2_1_c = zeros(length(t),3,length(nflex)); PHUB_3D = P2_1;
d_norma = zeros(length(t),1,length(nflex));
H_1 = zeros(length(t),3,length(nflex)); H_2 = H_1; H_net = zeros(length(t),3,length(nflex)); 
H_1_n = H_net; H_2_n = H_net; E_sys = H_net; H_bodies = zeros(length(t),3,length(nflex)); H_net_lin = H_net;
d_Ang = zeros(length(t),3,length(nflex));
domega = zeros(length(t),1,length(nflex)); domega_v = H_net; Mec = zeros(length(t),length(nflex));
CM = zeros(length(t),3);

%Evaluating all bodies
for j = 2:length(nflex)
    D2_1 = Bodies.B(j).Main_position;
    D1_2 = Bodies.B(j).Body_position;
    parfor i=1:length(t)
        %Pre-alocating variables
        K = []; 
        C321_K = [];
        vec_b1 = [];
        ArotB_I_BK = [];
        
        %Rotation matrix of Body 1 from the LVLH frame to the BRF
        rot = X(i, Bodies.B(1).pos_Theta)';
        if quat == 0
            phi = rot(1); theta = rot(2); psi = rot(3);
            C321_1 = angle2dcm(psi,theta,phi);
        else
            Sq = [0 -rot(3) rot(2); rot(3) 0 -rot(1); -rot(2) rot(1) 0];
            C321_1 = (rot(4)^2-rot(1:3)'*rot(1:3))*eye(3)+2*rot(1:3)*rot(1:3)'-2*rot(4)*Sq;
        end
        
        %Main body position vector in the ECI
        R1 = X(i,Bodies.B(1).pos_R); x = R1(1); y = R1(2); z = R1(3);
        rho = (x^2+y^2+z^2)^.5; delta = pi/2-acos(z/rho); lambda = atan2(y,x);
        %Rotation matrix from the ECI to the LVLH frame
        C32_B1 = angle2dcm(lambda,-delta-pi/2,0); 
        C32_B1 = angle2dcm(pi/2,0,0)*C32_B1;
        %Rotation matrix from ECI to BRF1
        ArotI_B_B1 = C321_1*C32_B1;
        %Rotation matrix from BRF1 to ECI
        ArotB_I_B1 = ArotI_B_B1';
        
        %Rotation matrix of Body 2 from the LVLH frame to the BRF2
        rot = X(i,Bodies.B(j).pos_Theta)';
        if quat == 0
            phi = rot(1); theta = rot(2); psi = rot(3);
            C321_2 = angle2dcm(psi,theta,phi);
        else
            Sq = [0 -rot(3) rot(2); rot(3) 0 -rot(1); -rot(2) rot(1) 0];
            C321_2 = (rot(4)^2-rot(1:3)'*rot(1:3))*eye(3)+2*rot(1:3)*rot(1:3)'-2*rot(4)*Sq;
        end
        
        %Second body position vector in the ECI
        R2 = X(i,Bodies.B(j).pos_R); x = R2(1); y = R2(2); z = R2(3);
        rho = (x^2+y^2+z^2)^.5; delta = pi/2-acos(z/rho); lambda = atan2(y,x);
        %Rotation matrix from the ECI to the LVLH frame
        C32_B2 = angle2dcm(lambda,-delta-pi/2,0); 
        C32_B2 = angle2dcm(pi/2,0,0)*C32_B2;
        %Rotation matrix from ECI to BRF2
        ArotI_B_B2 = C321_2*C32_B2;
        %Rotation matrix from BRF2 to ECI
        ArotB_I_B2 = ArotI_B_B2';
        
        %Second body position from Body 1 written in BRF1. The last term
        %substracts to generate an position error; thus all errors from the
        %booms positions can be displayes in the same graph.
        P2_1_c(i,:,j) = (ArotB_I_B1'*(X(i,Bodies.B(j).pos_R)'-X(i,Bodies.B(1).pos_R)') + ArotB_I_B1'*ArotB_I_B2*Bodies.B(j).Main_position - Bodies.B(j).Body_position)';
        d_norma(i,j) = norm(X(i,Bodies.B(1).pos_R)'-X(i,Bodies.B(j).pos_R)'-ArotB_I_B2*Bodies.B(j).Main_position);
        
        
        %PHUB_3D stores the relative position of body j regarding the HUB
        %in the LVLH ref system to be further used in the 3D animation
        if Bodies.System.Animation == 1
            PHUB_3D(i,:,j) = (C321_1'*((ArotB_I_B1'*(X(i,Bodies.B(j).pos_R)'-X(i,Bodies.B(1).pos_R)')+C321_1*C321_2'*Bodies.B(j).Main_position)))';
        end
        
        %Identify angle between bodies
        if Bodies.B(j).Constraint_type == 1
            if Bodies.B(j).Main_rotation==2
                vec_b1 = [1;0;0];
            else
                vec_b1 = [0;1;0];
            end
            if Bodies.B(j).Body_rotation == 2
                vec_b2 = [1;0;0];
            else
                vec_b2 = [0;1;0];
            end
            vec_b2_B1 = ArotI_B_B1*ArotB_I_B2*vec_b2;
            [dPsi, dTheta, dPhi] = dcm2angle(vrrotvec2mat(vrrotvec(vec_b1,vec_b2_B1)));
            d_Ang(i,:,j) = -[dPhi,dTheta,dPsi];
            domega(i,:,j) = (X(i,Bodies.B(j).pos_omega(Bodies.B(j).Main_rotation))-X(i,Bodies.B(1).pos_omega(Bodies.B(j).Body_rotation))*Bodies.B(j).Rot_Dir)';
            domega_v(i,:,j) = (C321_1*C321_2'*X(i,Bodies.B(j).pos_omega)'-X(i,Bodies.B(1).pos_omega)');
        else
            vec_b2_B1 = ArotI_B_B2*ArotB_I_BK*[0;1;0];
            [dPsi, dTheta, dPhi] = dcm2angle(vrrotvec2mat(vrrotvec([0; 1; 0],vec_b2_B1)));
            d_Ang(i,:,j) = -[dPhi, dTheta, dPsi];
            domega(i,:,j) = 0;
            domega_v(i,:,j) = (C321_K*C321_2'*X(i,Bodies.B(j).pos_omega)'-X(i,Bodies.B(K).pos_omega)');
        end
        
        r = Bodies.B(j).Main_position;
        if nflex(j) == 0
            Inertia = Bodies.B(j).Inertia;
        else
            Inertia = Bodies.B(j).Func.I_TT(X(i,Bodies.B(j).pos_qf)');
        end
        H_net(i,:,j) = ArotB_I_B2*(Inertia*X(i,Bodies.B(j).pos_omega)'); %In body 1 fixed body reference frame
        H_net_lin(i,:,j) = ((cross(X(i,Bodies.B(j).pos_R)',Bodies.B(j).mass*eye(3)*X(i,Bodies.B(j).pos_v)')))'; %In body 1 fixed body reference frame
        
        %Center of mass position at each time instant (must remain constant).
        if nflex(j) == 0
            CM(i,:) = CM(i,:) + Bodies.B(j).mass*(R2);
        else
            CM(i,:) = CM(i,:) + (Bodies.B(j).mass + Bodies.B(j).Tip_mass)*(R2 + (ArotB_I_B2*Bodies.B(j).Func.Center_of_mass(X(i,Bodies.B(j).pos_qf)'))');
        end
    end
end

%-------------------------------------------------------------------------------%
rot = X(end,Bodies.B(1).pos_Theta)';
if quat == 0
    phi = rot(1); theta = rot(2); psi = rot(3);
    C321_1 = angle2dcm(psi,theta,phi);
else
    Sq = [0 -rot(3) rot(2); rot(3) 0 -rot(1); -rot(2) rot(1) 0];
    C321_1 = (rot(4)^2-rot(1:3)'*rot(1:3))*eye(3)+2*rot(1:3)*rot(1:3)'-2*rot(4)*Sq;
end

R1 = X(end,Bodies.B(1).pos_R); x = R1(1); y = R1(2); z = R1(3);
rho = (x^2+y^2+z^2)^.5; delta = pi/2-acos(z/rho); lambda = atan2(y,x);
C32_B1 = angle2dcm(lambda,-delta-pi/2,0);
C32_B1 = angle2dcm(pi/2,0,0)*C32_B1;
ArotI_B_B1 = C321_1*C32_B1;
ArotB_I_B1 = ArotI_B_B1';
%-------------------------------------------------------------------------------%

for i = 2:length(nflex)
    disp(['Body ',num2str(i),' relative final position from central body - final: ']);
    disp(ArotB_I_B1'*(X(end,Bodies.B(i).pos_R)'-X(end,Bodies.B(1).pos_R)'));
    erro_x = abs(max(abs(P2_1_c(:,1,i)))/Bodies.B(i).Body_position(1)*100); %Maximum position errors
    erro_y = abs(max(abs(P2_1_c(:,2,i)))/Bodies.B(i).Body_position(2)*100); %Maximum position errors
    erro_z = abs(max(abs(P2_1_c(:,3,i)))/Bodies.B(i).Body_position(3)*100); %Maximum position errors
    disp(['Position error in x: ',num2str(erro_x),'%.']);
    disp(['Position error in y: ',num2str(erro_y),'%.']);
    disp(['Position error in z: ',num2str(erro_z),'%.']);
    E_abs = max(abs(d_norma(:,1,i)-d_norma(1,1,i)))/d_norma(1,1,i)*100;
    disp(['Absolute error: ',num2str(E_abs),'%.']);
    Bodies.SimOutput.B(i).erro = [erro_x erro_y erro_z];
    Bodies.SimOutput.B(i).E_abs = E_abs;
end


for j=2:length(nflex)
    figure(nf)
    subplot(311); plot(t,P2_1_c(:,1,j)); xlabel(' t (s) '); ylabel(' x (m) '); grid minor;
    subplot(312); plot(t,P2_1_c(:,2,j)); xlabel(' t (s) '); ylabel(' y (m) '); grid minor;
    subplot(313); plot(t,P2_1_c(:,3,j)); xlabel(' t (s) '); ylabel(' z (m) '); grid minor;
    sgtitle(['Position error ', Bodies.B(j).Name ,' (Body ',num2str(j),') and Body 1 - BRF^1']);
    nf=nf + 1;
end

if length(nflex)>1
    clear lgC;
    figure(nf)
    for j=2:length(nflex)
        plot(t,d_norma(:,1,j)); xlabel(' t (s) '); ylabel(' Distance (m) ');
        title('Bodies distance - norm'); hold on;
        lgC{j-1}=num2str(j,'Body %-.0f');
    end
    legend(lgC); grid minor;
    nf = nf + 1;
end


for j = 2:length(nflex)
    figure(nf)
    subplot(321); plot(t,d_Ang(:,1,j)*180/pi); xlabel(' t (s) '); ylabel(' d\phi (degrees) '); grid minor;
    subplot(322); plot(t,domega_v(:,1,j)*180/pi); xlabel(' t (s) '); ylabel(' p (degrees/s) '); grid minor;
    subplot(323); plot(t,d_Ang(:,2,j)*180/pi); xlabel(' t (s) '); ylabel(' d\theta (degrees) '); grid minor;
    subplot(324); plot(t,domega_v(:,2,j)*180/pi); xlabel(' t (s) '); ylabel(' q (degrees/s) '); grid minor;
    subplot(325); plot(t,d_Ang(:,3,j)*180/pi); xlabel(' t (s) '); ylabel(' d\psi (degrees) '); grid minor;
    subplot(326); plot(t,domega_v(:,3,j)*180/pi); xlabel(' t (s) '); ylabel(' r (degrees/s) '); grid minor;
    sgtitle(['Relative angular position and velocity between', Bodies.B(j).Name, ' (Body ',num2str(j),') and 1 - BRF^1']);
    nf = nf + 1;
end


% Code section used to evaluate the system's angular momentum
% and CM position.
H_netN = zeros(length(t),1); MecN=H_netN; H_netN_lin=H_netN;
Inertia = Bodies.B(1).Inertia;
pos_omega = Bodies.B(1).pos_omega;
pos_v = Bodies.B(1).pos_v;
pos_R = Bodies.B(1).pos_R;
for i = 1:length(t)
   rot = X(i,Bodies.B(1).pos_Theta)';
   if quat == 0 
       phi = rot(1); theta = rot(2); psi = rot(3);
       C321_1 = angle2dcm(psi,theta,phi);
   else
       Sq = [0 -rot(3) rot(2); rot(3) 0 -rot(1); -rot(2) rot(1) 0];
       C321_1 = (rot(4)^2-rot(1:3)'*rot(1:3))*eye(3)+2*rot(1:3)*rot(1:3)'-2*rot(4)*Sq;
   end
   R1 = X(i,Bodies.B(1).pos_R); x = R1(1); y = R1(2); z = R1(3);
   rho = (x^2+y^2+z^2)^.5; delta = pi/2-acos(z/rho); lambda = atan2(y,x);
   C32_B1 = angle2dcm(lambda,-delta-pi/2,0);
   C32_B1 = angle2dcm(pi/2,0,0)*C32_B1;
   ArotI_B_B1 = C321_1*C32_B1; 
   R1_B1 = ArotI_B_B1*R1';
   Inertia = Bodies.B(1).Inertia;
    
   H_net(i,:,1) = ArotI_B_B1'*Inertia*X(i,pos_omega)';
   H_net_lin(i,:,1) = ((cross(X(i,pos_R)',Bodies.B(1).mass*eye(3)*X(i,pos_v)')))';
   
   temp = zeros(3,1);
   temp_lin = temp;
   for j = 1:length(nflex)
       temp = temp+H_net(i,:,j)';
       temp_lin = temp_lin+H_net_lin(i,:,j)';
   end
   H_netN(i) = norm(temp);
   H_netN_lin(i) = norm(temp_lin+temp);
   
   CM(i,:) = CM(i,:) + Bodies.B(1).mass*R1;
   CM_n(i) = norm(CM(i,:));
end

% For to evaluate the system's total mass in order to obtain the system's
% center o mass (CM).
M_T = 0;
for i = 1:length(nflex)
    M_T = M_T + Bodies.B(i).mass;
end

CM = CM/M_T;


% Center of mass should be constant for a circular orbit without applied
% forces in the system
figure(nf)
subplot(211); plot3(CM(:,1), CM(:,2), CM(:,3)); xlabel(' x (m) '); ylabel(' y (m) '); ylabel(' z (m) '); grid minor;
subplot(212); plot(t, CM_n - CM_n(1)); xlabel(' t (s) '); ylabel(' Norm CG error (m) '); grid minor;
sgtitle('Systems Center of mass - ECI'); 
nf = nf + 1;


figure(nf)
plot(t,H_netN); xlabel(' t (s) '); ylabel(' kgm^2/s'); grid minor;
sgtitle('Systems Angular Momentum - BRF1'); legend('H_{net}');
nf = nf + 1;

figure(nf)
plot(t,H_netN_lin); xlabel(' t (s) '); ylabel(' kgm^2/s'); grid minor;
sgtitle('Systems Total Angular Momentum - ECI'); legend('H_{net} Total');
nf = nf + 1;


%Function to analyze the vibration of the booms
if Bodies.System.FFT == 1
    if Bodies.System.LFFT < max(t)
        if Bodies.System.LFFT > 0
            if max(t)-Bodies.System.LFFT > 2
                for j = 2:length(nflex)
                    figure(nf)
                    subplot(321); plot(t,d_Ang(:,1,j)*180/pi); xlabel(' t (s) '); ylabel(' d\phi (degrees) '); grid minor; xlim([Bodies.System.LFFT Bodies.System.LFFT+1]);
                    subplot(322); plot(t,domega_v(:,1,j)*180/pi); xlabel(' t (s) '); ylabel(' p (degrees/s) '); grid minor; xlim([Bodies.System.LFFT Bodies.System.LFFT+1]);
                    subplot(323); plot(t,d_Ang(:,2,j)*180/pi); xlabel(' t (s) '); ylabel(' d\theta (degrees) '); grid minor; xlim([Bodies.System.LFFT Bodies.System.LFFT+1]);
                    subplot(324); plot(t,domega_v(:,2,j)*180/pi); xlabel(' t (s) '); ylabel(' q (degrees/s) '); grid minor; xlim([Bodies.System.LFFT Bodies.System.LFFT+1]);
                    subplot(325); plot(t,d_Ang(:,3,j)*180/pi); xlabel(' t (s) '); ylabel(' d\psi (degrees) '); grid minor; xlim([Bodies.System.LFFT Bodies.System.LFFT+1]);
                    subplot(326); plot(t,domega_v(:,3,j)*180/pi); xlabel(' t (s) '); ylabel(' r (degrees/s) '); grid minor; xlim([Bodies.System.LFFT Bodies.System.LFFT+1]);
                    sgtitle(['Relative angular position and velocity between', Bodies.B(j).Name, ' (Body ',num2str(j),') and 1 - BRF^1: FFT INITIAL TIME INTERVAL']);
                    nf = nf + 1;
                end
            end
            for j = 2:length(nflex)
                figure(nf)
                subplot(321); plot(t,d_Ang(:,1,j)*180/pi); xlabel(' t (s) '); ylabel(' d\phi (degrees) '); grid minor; xlim([Bodies.System.LFFT max(t)]);
                subplot(322); plot(t,domega_v(:,1,j)*180/pi); xlabel(' t (s) '); ylabel(' p (degrees/s) '); grid minor; xlim([Bodies.System.LFFT max(t)]);
                subplot(323); plot(t,d_Ang(:,2,j)*180/pi); xlabel(' t (s) '); ylabel(' d\theta (degrees) '); grid minor; xlim([Bodies.System.LFFT max(t)]);
                subplot(324); plot(t,domega_v(:,2,j)*180/pi); xlabel(' t (s) '); ylabel(' q (degrees/s) '); grid minor; xlim([Bodies.System.LFFT max(t)]);
                subplot(325); plot(t,d_Ang(:,3,j)*180/pi); xlabel(' t (s) '); ylabel(' d\psi (degrees) '); grid minor; xlim([Bodies.System.LFFT max(t)]);
                subplot(326); plot(t,domega_v(:,3,j)*180/pi); xlabel(' t (s) '); ylabel(' r (degrees/s) '); grid minor; xlim([Bodies.System.LFFT max(t)]);
                sgtitle(['Relative angular position and velocity between', Bodies.B(j).Name, ' (Body ',num2str(j),') and 1 - BRF^1: FFT TIME INTERVAL']);
                nf = nf + 1;
            end 
        end
        fft_study
    else
        disp('Error: Initial FFT time is greater than simulation time');
    end
end

global d_Ang_ord Ang_ord PHUB_3D_ord
if Bodies.System.Animation == 1
    FR = Bodies.System.FR; %Sampling frequency
    sample = 1/FR;
    d_Ang_ord = zeros(FR*max(t),3,length(nflex)-1);
    PHUB_3D_ord = d_Ang_ord;
    for i = 2:length(nflex)
        teste = d_Ang(:,:,i);
        testeHUB = PHUB_3D(:,:,i);
        
        %Starts the GPU
        gpuDevice(1);
        %Allocates vectors in GPU
        t_GPU = gpuArray(t);
        t_vector = 0:sample:(max(t)-1*sample); t_vector = gpuArray(t_vector);
        Vec_ORD = gpuArray(sortrows([t, teste, testeHUB]));
        
        %Interpolate the results
        d_Ang_ord(:,1,i) = gather(interp1(t_GPU,Vec_ORD(:,2),t_vector, 'spline'));
        d_Ang_ord(:,2,i) = gather(interp1(t_GPU,Vec_ORD(:,3),t_vector, 'spline'));
        d_Ang_ord(:,3,i) = gather(interp1(t_GPU,Vec_ORD(:,4),t_vector, 'spline'));
        PHUB_3D_ord(:,1,i) = gather(interp1(t_GPU,Vec_ORD(:,5),t_vector, 'spline'));
        PHUB_3D_ord(:,2,i) = gather(interp1(t_GPU,Vec_ORD(:,6),t_vector, 'spline'));
        PHUB_3D_ord(:,3,i) = gather(interp1(t_GPU,Vec_ORD(:,7),t_vector, 'spline'));
        
    end
    Ang_ord = zeros(FR*max(t)-1,3,length(nflex));
    for i = 1:length(nflex)
        vec = X(:,Bodies.B(i).pos_Theta);
        %Starts the GPU
        gpuDevice(1);
        %Allocates vectors in GPU
        t_GPU = gpuArray(t);
        t_vector = 0:sample:(max(t)-1*sample); t_vector = gpuArray(t_vector);
        Vec_ORD = gpuArray(sortrows([t, vec]));
        rot1 = gather(interp1(t_GPU,Vec_ORD(:,2),t_vector, 'spline'))';
        rot2 = gather(interp1(t_GPU,Vec_ORD(:,3),t_vector, 'spline'))';
        rot3 = gather(interp1(t_GPU,Vec_ORD(:,4),t_vector, 'spline'))';
        rot4 = gather(interp1(t_GPU,Vec_ORD(:,5),t_vector, 'spline'))';
        rot_full = [rot1, rot2, rot3, rot4];
        for j = 1:FR*max(t)

            rot = rot_full(j,:)';
            Sq = [0 -rot(3) rot(2); rot(3) 0 -rot(1); -rot(2) rot(1) 0];
            C321 = (rot(4)^2-rot(1:3)'*rot(1:3))*eye(3)+2*rot(1:3)*rot(1:3)'-2*rot(4)*Sq;
            [Psi, Theta, Phi] = dcm2angle(C321);
            
            
            Ang_ord(j,1,i) = Phi;
            Ang_ord(j,2,i) = Theta;
            Ang_ord(j,3,i) = Psi;
        end
    end
    %Load World data
    world = vrworld('Modelo_3D');
    open(world);
    % view(world);
    sim('Animation_Simulink',[0 max(t)-5/FR]);
end