function in=initial(Y)
%Função para verificar as condições iniciais

global X0 quat nflex Bodies

%Inicialmente assume-se que X=X0 determinado pelo usuário
X=X0;

%Aloca os vetores de Y no vetor inicial X
X(7:length(Y)+6)=Y;

%Obtém a derivada dos estados:
qpp=dinam_sat_flex(0,X);

%Equações de restrição para o problema - devem ser respeitadas na condição
%inicial

%Alocando os estados
% rot=X(end-3*length(nflex)+1:end);
rot=X(Bodies.B(1).pos_Theta(1):end);
omega=[];
for i=1:length(nflex)
    el=6*i+sum(nflex(1:i))-(6+nflex(i))+1;
    omega(:,i)=X(el+3:el+6-1);
    aux=[0,cumsum(nflex)];
    el=length(nflex)*6+sum(nflex)+3*(i-1)+aux(i)+1;
    R(:,i)=X(el:el+3-1);
end

%Matrizes de rotação para o problema:
C321=zeros(3,3,length(nflex)); 
for i=1:length(nflex)
    if quat==0 %Caso utilize angulos de Euler
        phi=rot(1+3*(i-1)); theta=rot(2+3*(i-1)); psi=rot(3+3*(i-1));
        C321(:,:,i)=angle2dcm(psi,theta,phi); %Mudei aqui %Conferir a matriz
    else %Caso uilize quaternions
        roti=rot(1+4*(i-1):4+4*(i-1));
        Sq=[0 -roti(3) roti(2); roti(3) 0 -roti(1); -roti(2) roti(1) 0];
        C321(:,:,i)=(roti(4)^2-roti(1:3)'*roti(1:3))*eye(3)+2*roti(1:3)*roti(1:3)'-2*roti(4)*Sq;
    end
end

C321=zeros(3,3,length(nflex)); C32=zeros(3,3,length(nflex)); ArotI_B=zeros(3,3,length(nflex)); ArotB_I=zeros(3,3,length(nflex));
for i=1:length(nflex)
    if quat==0 %Caso utilize angulos de Euler
        phi=rot(1+3*(i-1)); theta=rot(2+3*(i-1)); psi=rot(3+3*(i-1));
        C321(:,:,i)=angle2dcm(psi,theta,phi); %Mudei aqui %Conferir a matriz
    else %Caso uilize quaternions
        roti=rot(1+4*(i-1):4+4*(i-1));
        Sq=[0 -roti(3) roti(2); roti(3) 0 -roti(1); -roti(2) roti(1) 0];
        C321(:,:,i)=(roti(4)^2-roti(1:3)'*roti(1:3))*eye(3)+2*roti(1:3)*roti(1:3)'-2*roti(4)*Sq;
    end
    %Dados de posição
    x=R(1,i); y=R(2,i); z=R(3,i); rho=(x^2+y^2+z^2)^.5;
    %Angulos em coordenadas esféricas do sistema LVLH em relação ao ECI
    delta=pi/2-acos(z/rho); lambda=atan2(y,x);
    C32_B=angle2dcm(lambda,-delta-pi/2,0);  %Conferir a matriz
    C32(:,:,i)=angle2dcm(pi/2,0,0)*C32_B;
    ArotI_B(:,:,i)=C321(:,:,i)*C32(:,:,i);
    ArotB_I(:,:,i)=ArotI_B(:,:,i)';
end

Cq_f=[]; Qc_f=[];
Q1=quaternion(ArotB_I(:,:,1)');
for i=2:length(nflex)
    Qi=quaternion(ArotB_I(:,:,i)');
    Cq=zeros(6,length(nflex)*6+sum(nflex));
    [Gb_1,k_1,Cqf_1]=constraints_trl_quat(omega(:,1),Q1,Bodies.B(i).Body_position);
    [Gb_i,k_i,Cqf_i]=constraints_trl_quat(omega(:,i),Qi,Bodies.B(i).Main_position);
    [Cq_1,Cq_i,Cqf]=constraints_rot_nova_quat(omega(:,1),omega(:,i),Q1,Qi,i,Bodies);
    
    %Número dos estados do corpo   
    tr_i=6*(i-1)+cumsum(nflex(i-1))+1; %Posição dos estados de translação
    rot_i=tr_i+3;
    
    %Restrição de translação
    Cq(1:3,1:3)=-eye(3);   %Sempre nessas posições, pois é o primeiro corpo
    Cq(1:3,4:6)=-k_1*Gb_1; %Sempre nessas posições, pois é o primeiro corpo
    Cq(1:3,tr_i:tr_i+2)=eye(3);
    Cq(1:3,rot_i:rot_i+2)=k_i*Gb_i;
    %Restrição de rotação
    Cq(4:6,4:6)=Cq_1; %Sempre nessas posições, pois é o primeiro corpo
    Cq(4:6,rot_i:rot_i+2)=Cq_i;
    
    %Forças das restrições de translação
    Qc(1:3,1)=-Cqf_i-(-Cqf_1); 
    %Forças das restrições de rotação
    Qc(4:6,1)=-Cqf;
    
    Cq_f=[Cq_f;Cq]; Qc_f=[Qc_f;Qc];
end

%O que deve-se manter zero durante todos os instantes:
%Talvez adicionar a derivada primeira da equação de restrição se não der
%certo
%First Lagrange relation
qpp_dim=6*length(nflex)+sum(nflex);
l_1=Qc_f-Cq_f*qpp(1:qpp_dim);
l_2=Cq_f*X(1:qpp_dim);

in=[l_1
    l_2];