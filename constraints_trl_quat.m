function [G, k, Cqf] = constraints_trl_quat(omega, ang_i, D1_2)
%//////////////////////////////////////////////////////////////////%
%Function used to generate translation constraints (constraints three
%degrees of freedom)
%See: Shabana, A. A. Computational Dynamics. ISBN: 9780470686850. DOI: 10.1002/9780470686850
%INPUTS:                              
%omega: Angular velocity in the BRF (i) from the ECI
%ang_i: Attitude quaternions from the ECI to BRF (i)
%D1_2: Attached body (j) distance in the BRF (i)
%OUTPUTS:
%G: Matriz relating angular velocity vector and quaternions time derivative
%k: Matrix of partial derivatives
%Cqf: 'Force' vector associated with the constraints
%//////////////////////////////////////////////////////////////////%
u01 = D1_2(1); u02 = D1_2(2); u03 = D1_2(3);

%Quaternion components
Q1 = ang_i(1); Q2 = ang_i(2); Q3 = ang_i(3); Q4 = ang_i(4);

%Angular velocity in the BRF from the ECI
p = omega(1); q = omega(2); r = omega(3);

k = [2.*Q1.*u01+2.*Q2.*u02+2.*Q3.*u03,(-2).*Q2.*u01+2.*Q1.*u02+ ...
  2.*Q4.*u03,(-2).*Q3.*u01+(-2).*Q4.*u02+2.*Q1.*u03,2.*Q4.* ...
  u01+(-2).*Q3.*u02+2.*Q2.*u03;2.*Q2.*u01+(-2).*Q1.*u02+(-2).* ...
  Q4.*u03,2.*Q1.*u01+2.*Q2.*u02+2.*Q3.*u03,2.*Q4.*u01+(-2).* ...
  Q3.*u02+2.*Q2.*u03,2.*Q3.*u01+2.*Q4.*u02+(-2).*Q1.*u03;2.* ...
  Q3.*u01+2.*Q4.*u02+(-2).*Q1.*u03,(-2).*Q4.*u01+2.*Q3.*u02+( ...
  -2).*Q2.*u03,2.*Q1.*u01+2.*Q2.*u02+2.*Q3.*u03,(-2).*Q2.*u01+ ...
  2.*Q1.*u02+2.*Q4.*u03];

G = [(1/2).*Q4,(-1/2).*Q3,(1/2).*Q2;(1/2).*Q3,(1/2).*Q4,(-1/2).* ...
  Q1;(-1/2).*Q2,(1/2).*Q1,(1/2).*Q4;(-1/2).*Q1,(-1/2).*Q2,( ...
  -1/2).*Q3];

Cq_qp_q = [2.*((-1).*q.*Q3.*u01+Q2.*r.*u01+p.*Q3.*u02+(-1).*Q1.*r.* ...
  u02+q.*Q1.*u03+(-1).*p.*Q2.*u03),2.*(Q1.*r.*u01+Q2.*r.*u02+ ...
  Q4.*((-1).*q.*u01+p.*u02)+(-1).*p.*Q1.*u03+(-1).*q.*Q2.*u03) ...
  ,2.*((-1).*Q4.*r.*u01+p.*Q1.*u02+Q3.*r.*u02+p.*Q4.*u03+(-1) ...
  .*q.*(Q1.*u01+Q3.*u03)),2.*((-1).*q.*Q2.*u01+(-1).*Q3.*r.* ...
  u01+p.*Q2.*u02+(-1).*Q4.*r.*u02+p.*Q3.*u03+q.*Q4.*u03);2.*( ...
  q.*Q4.*u01+(-1).*p.*Q4.*u02+(-1).*r.*(Q1.*u01+Q2.*u02)+p.* ...
  Q1.*u03+q.*Q2.*u03),2.*((-1).*q.*Q3.*u01+Q2.*r.*u01+p.*Q3.* ...
  u02+(-1).*Q1.*r.*u02+q.*Q1.*u03+(-1).*p.*Q2.*u03),2.*((-1).* ...
  q.*Q2.*u01+(-1).*Q3.*r.*u01+p.*Q2.*u02+(-1).*Q4.*r.*u02+p.* ...
  Q3.*u03+q.*Q4.*u03),2.*((-1).*(p.*Q1+Q3.*r).*u02+Q4.*(r.* ...
  u01+(-1).*p.*u03)+q.*(Q1.*u01+Q3.*u03));2.*((-1).*(p.*Q1+ ...
  Q3.*r).*u02+Q4.*(r.*u01+(-1).*p.*u03)+q.*(Q1.*u01+Q3.*u03)), ...
  2.*(q.*Q2.*u01+Q3.*r.*u01+(-1).*p.*Q2.*u02+(-1).*p.*Q3.*u03+ ...
  Q4.*(r.*u02+(-1).*q.*u03)),2.*((-1).*q.*Q3.*u01+Q2.*r.*u01+ ...
  p.*Q3.*u02+(-1).*Q1.*r.*u02+q.*Q1.*u03+(-1).*p.*Q2.*u03),2.* ...
  (Q1.*r.*u01+Q2.*r.*u02+Q4.*((-1).*q.*u01+p.*u02)+(-1).*p.* ...
  Q1.*u03+(-1).*q.*Q2.*u03)];

Cqf = Cq_qp_q*G*omega;