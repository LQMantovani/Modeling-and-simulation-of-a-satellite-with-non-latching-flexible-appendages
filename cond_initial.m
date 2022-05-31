function [R, V, omega, qf, qfp, Bodies]=cond_initial(i, Bodies)
%Function to generate the initial conditions for each body in the system
%i: Body number
%Bodies: Data structure with bodies information


%Orbital parameters in spheric coordinate system
%R_0: Orbit radius
%Omega(1): lambda (longitude)
%Omega(2): delta (latitude)
%Theta: Angular position of body i, in Euler angles (phi,theta,psi)
%nflex: Number of elastic modes
%D1_2: Distance from body i to main body
R_0 = Bodies.System.Radius;
lambda = Bodies.System.Lambda;    %Central body Longitude
delta = Bodies.System.Delta;      %Central body Latitudnoe
muT = Bodies.System.muT;
Theta_0 = Bodies.B(i).Angular_pos;
nflex = Bodies.B(i).nflex;
Di_j = Bodies.B(i).Main_position;
Dj_i = Bodies.B(i).Body_position;


%Initial position
%Rotation matrix - from the Earth centered inertial frame (ECI) to the LVLH
C32 = angle2dcm(lambda,-delta-pi/2,0);
C32 = angle2dcm(pi/2,0,0)*C32;
%Rotation matrix from the LVLH to the body fixed reference frame
%Booms angular position are given in relation to the HUB in
%bodies_data.txt. This is the reason to use this if.
if i>1
    C321 = angle2dcm(Theta_0(3),Theta_0(2),Theta_0(1))*Bodies.B(1).C321; %Rotation matrix from the LVLH to the body fixed frame
else
    C321 = angle2dcm(Theta_0(3),Theta_0(2),Theta_0(1)); %Rotation matrix deom the LVLH to the body fixed frame
end
ArotI_B = C321*C32; %Rotation matrix from the inertia system to the body fixed frame
ArotB_I = ArotI_B'; %Rotation matrix from the body fixed frame to the ECI frame


Bodies.B(i).C32 = C32;
Bodies.B(i).C321 = C321;
Bodies.B(i).ArotB_I = ArotB_I;

%Position in the ECI frame
%Position vector in the inertial system:
if any([i == 1, Bodies.B(i).Constraint_type == 1])
    R_I = C32'*[0;0;-R_0]+Bodies.B(1).C32'*(Bodies.B(1).C321'*Dj_i+C321'*-Di_j);
else
    j = Bodies.B(i).Attached_body;
    vec_pos = zeros(3,1); vec_pos(Bodies.B(j).Axis_Length) = Bodies.B(j).Length;
    R_I = Bodies.B(j).R_I + (Bodies.B(j).ArotB_I*vec_pos);
end

Bodies.B(i).R_I = R_I;

if i>1
    %Recalculating the inertial position angles lambda and delta now that the
    %position in the inertial system is known
    rho = norm(R_I);
    delta = pi/2-acos(R_I(3)/rho); lambda=atan2(R_I(2),R_I(1));
    C32 = angle2dcm(lambda,-delta-pi/2,0);
    C32 = angle2dcm(pi/2,0,0)*C32;
    ArotI_B = C321*C32; 
    ArotB_I = ArotI_B'; 
    Bodies.B(i).C32 = C32;
    Bodies.B(i).ArotB_I = ArotB_I;
end


n = -sqrt(muT/R_0^3);  %Orbit frequency (circular orbit) in the LVLH system
omega_LVLH = [0;-n;0]; %Angular velocity in the LVLH frame
omega_0_B = C321*(omega_LVLH+Bodies.B(1).C321'*Bodies.B(1).Angular_vel); %Angular velocity in the body fixed frame

%Velocity in the ECI frame
if any([i == 1, Bodies.B(i).Constraint_type == 1])
    V0 = cross(C32'*omega_LVLH,R_I)+C32'*cross(Bodies.B(1).C321'*Bodies.B(1).Angular_vel,Bodies.B(1).C321'*Dj_i+C321'*-Di_j);
else
    V0 = cross(C32'*omega_LVLH,R_I)+C32'*cross(Bodies.B(1).C321'*Bodies.B(1).Angular_vel,Bodies.B(1).C321'*Dj_i+Bodies.B(j).C321'*(-Di_j + vec_pos));
end

%Initial conditions for body i
V = V0;               %Velocity in the ECI frame
omega = omega_0_B;     %Angular velocity in the body i fixed frame
qfp = zeros(nflex,1); %Elastic modes velocity in the body i fixed frame
R = R_I;              %Position in the ECI frame
qf = zeros(nflex,1);  %Elastic modes in the body i fixed frame 