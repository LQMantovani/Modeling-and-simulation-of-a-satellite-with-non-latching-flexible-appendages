%//////////////////////////////////////////////////////////////////%
%MAIN FILE
%//////////////////////////////////////////////////////////////////%

clear all; close all; clc;

global var

%Read the input data:
filename = 'bodies_data.txt'; %Used to generate the Bodies_data.mat para o caso do Integrated SPORT - Parallel
data_read

%Determine the Generalized matrices for each body when necessary
Allocate_matrices

%System information   
raioT=6378.135e3;               %Equatorial radius in meters               
Gu=6.67408e-11;                 %Gravitacional constant
mT=5.973332e24;                 %Earth's mass
fr=pi/180;                      %Convertion from degree to radian
Bodies.System.Radius = raioT + Bodies.System.Altitude;  %Orbit radius
Bodies.System.muT = Gu*mT;                              %Earth's mu
Bodies.System.nflex = nflex;

%Used to verify if the correct input is given
if Bodies.System.quat ~= 1, Bodies.System.quat = 0; end %If quat==1, quaternions are used
quat = Bodies.System.quat;

%User input with desired simulation time
Bodies.System.tsim=input('Simulation time (s): ');


%INITIAL CONDITIONS - @cond_initial
qp=[]; q=[]; Theta=[];
for i = 1:size(Bodies.B,2)
    [R_0, V_0, Omega_0, qf_0, qfp_0, Bodies] = cond_initial(i, Bodies);
    qp = [qp;V_0;Omega_0;qfp_0];
    q = [q; R_0; qf_0];
    if Bodies.System.quat==0
        Theta = [Theta; Bodies.B(i).Angular_pos];
    else
        qt = quaternion(Bodies.B(i).C321);
        Theta = [Theta; qt(1); qt(2); qt(3); qt(4)];
    end
end


%Defines the total number of constraint in the system based on the
%constraint type
n_constraints = 0;
for i = 2:length(nflex)
   if Bodies.B(i).Constraint_type == 1
       if Bodies.System.Lock
            n_constraints = n_constraints + 6;
       else
           n_constraints = n_constraints + 5;
       end
   else
       n_constraints = n_constraints + 6;
   end
end

Bodies.System.Constraint_number = n_constraints;

%Update states and states number
X0 = [qp;q;Theta];
n_estados = length(X0);


%To display the initial condition of each body
disp('----------------------------------------------------------------------');
for i = 1:length(nflex)
   disp(['Body ',num2str(i),' velocity in the ECI (m/s): ']);
   disp(X0(Bodies.B(i).pos_v));
   disp(['BOdy ',num2str(i),' angular velocity in the body reference frame with respected to the ECI (degrees/s):'])
   disp(X0(Bodies.B(i).pos_omega)*180/pi);
   if i>1
       disp(['Body ',num2str(i),' relative position from the main body:']);
       disp(Bodies.B(1).ArotB_I'*(X0(Bodies.B(i).pos_R) - X0(Bodies.B(1).pos_R)));
   end
end


%TO PERFORM THE INTEGRATION
%Setting integrator parameters as defined in the input file
opt_fsolve = odeset('RelTol', Bodies.System.RelTol, 'AbsTol', Bodies.System.AbsTol, 'MaxStep', Bodies.System.MaxStep);

var = .05; %Variable used to verify percentage of simulation
tic %To time the system - begin
[t, X] = ode15s(@dinam_sat_flex, [0 Bodies.System.tsim], X0, opt_fsolve, Bodies);
disp('Simulation 100% complete');
toc %To time the system - end

%TO GENERATE PLOTS
graphics