function X_flex=dinam_sat_flex(t, X, Bodies)
%//////////////////////////////////////////////////////////////////%
%Function used to organize variables and angular position equations 
%for each body
%Inputs:
%t: time
%X: State vector
%Bodies: Data structure with bodies data
%Outputs:
%X_flex: %States time derivative
%//////////////////////////////////////////////////////////////////%


global var
muT = Bodies.System.muT;
nflex = Bodies.System.nflex;
tsim = Bodies.System.tsim;

if t/tsim >= var
    disp(['Simulation ',num2str(floor(t/tsim*100)),'% complete']);
    var = var+0.05;
end
    

%Alocate the variables considering the number of bodies and number of
%flexible modes of each
X_flex_p = []; R_p = zeros(3,length(nflex)); R = R_p; omega = R_p;
for i = 1:length(nflex)
    el = 6*i+sum(nflex(1:i))-(6+nflex(i))+1;
    R_p(:,i) = X(el:el+3-1);
    omega(:,i) = X(el+3:el+6-1);
    if nflex(i) ~= 0
        B(i).Qfp = X(el+6:el+6+nflex(i)-1);
        X_flex_p = [X_flex_p; R_p(:,i); B(i).Qfp];
    else
        X_flex_p = [X_flex_p; R_p(:,i)];
    end
end

%Alocate variables considering the number of bodies and number of flexibles
%modes of each - this variables have already been integrated once
for i = 1:length(nflex)
    aux = [0, cumsum(nflex)];
    el = length(nflex)*6+sum(nflex)+3*(i-1)+aux(i)+1;
    R(:,i) = X(el:el+3-1);
    if nflex(i) ~= 0  
        B(i).Qf = X(el+3:el+3+nflex(i)-1); 
    else
       B(i).Qf = []; 
    end
end

quat = Bodies.System.quat;
rot = X(end-(3+quat)*length(nflex)+1:end);

%Normalizing the quaternion vectors
for i = 1:length(nflex)
    rot(1+4*(i-1):4+4*(i-1))=rot(1+4*(i-1):4+4*(i-1))/norm(rot(1+4*(i-1):4+4*(i-1)));
end

U = control(t, X, Bodies);

%Obtains both the generalized mass matrix and the generalized forces vector of the system:
[Ml, Ql] = eq_sistem(omega, rot, R, B, Bodies, R_p, U);


if Bodies.System.Embedded == false
    Xl = Ml\Ql;
    X_flex_pp = Xl(1:size(Xl,1) - Bodies.System.Constraint_number);
else
    X_flex_pp = Ml;
end



%Solving the angular positions derivatives for each body:
rotp = [];
for i = 1:length(nflex)
    if quat == 0
        n = -sqrt(muT/ norm(R(:,i))^3);
        p = omega(1,i); q = omega(2,i); r = omega(3,i);
        phi = rot(1+3*(i-1)); theta = rot(2+3*(i-1)); psi = rot(3+3*(i-1));
        if theta == pi/2, theta = pi/2-pi/2*.0001; end %To avoid the singularity ad theta = pi/2
        %Obtains the Euler angles time derivative
        phip = p+n*sec(theta)*sin(psi)+(cos(phi)*r+q*sin(phi))*tan(theta);
        thetap = n*cos(psi)+cos(phi)*q-r*sin(phi);
        psip = sec(theta)*(cos(phi)*r+q*sin(phi))+n*sin(psi)*tan(theta);
        rotp = [rotp; phip; thetap; psip];
    else
        n = -sqrt(muT/ norm(R(:,i))^3);
        roti = rot(1+4*(i-1):4+4*(i-1));
        Sq = [0 -roti(3) roti(2); roti(3) 0 -roti(1); -roti(2) roti(1) 0];
        C321 = (roti(4)^2-roti(1:3)'*roti(1:3))*eye(3)+2*roti(1:3)*roti(1:3)'-2*roti(4)*Sq;
        omega(:,i) = omega(:,i)-C321*[0;-n;0];
        %Generating the quaternions derivative. See: Tewari, A. Atmospheric
        %and Space Flight Dynamics. DOI: 10.1007/978-0-8176-4438-3
        p = omega(1,i); q = omega(2,i); r = omega(3,i);
        OMEGA = [0 r -q p; -r 0 p q; q -p 0 r; -p -q -r 0];
        rotp = [rotp; .5*OMEGA*roti];
    end
end

X_flex = [X_flex_pp; X_flex_p; rotp];
