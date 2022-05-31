function Xq = plot_degree(Xd)
%Xq = plot_degree(Xd)
% Function to obtain Euler Angles from attitude quaternions
% INPUT:
% Xd: Quaternion
% OUTPUT:
% Euler Angle

% Rotation matrix obtained from quaternion.
% See: Tewari, A. Atmospheric and Space Flight Dynamics. DOI: 10.1007/978-0-8176-4438-3

Xq = zeros(size(Xd,1), 3, 1);

for j = 1:size(Xd,1)
    rot = Xd(j,:)';
    Sq = [0 -rot(3) rot(2); rot(3) 0 -rot(1); -rot(2) rot(1) 0];
    C321 = (rot(4)^2-rot(1:3)'*rot(1:3))*eye(3)+2*rot(1:3)*rot(1:3)'-2*rot(4)*Sq;
    [Psi, Theta, Phi] = dcm2angle(C321);
    Xq(j,:,1) = [Phi; Theta; Psi]';
end