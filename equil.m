function eq=equil(X)
%Função para obter o equilíbrio
global Xe Ue nflex
%Velocidade angular
p=X(1); q=X(2); r=X(3);
X_t=Xe;

X_t(4)=p; X_t(10)=p;
X_t(5)=q; X_t(11)=q;
X_t(6)=r; X_t(12)=r;
%Angulos de Euler
X_t(end-2:end)=X(4:6);
X_t(end-5:end-3)=X_t(end-2:end);
%Coordenadas deformadas - velocidades

Xp=dinam_sat_flex(0,X_t);

%omega3p rotp qf1p qf2p qf1p qf2p

eq=[Xp(4:6);Xp(10:12);Xp(end-5:end)];