function [sys,x0,str,ts,simStateCompliance] = S_Function_Animation(t,x,u,flag)
global Bodies Ang_ord PHUB_3D_ord
switch flag
  case 0
       sizes = simsizes;

sizes.NumContStates  = 0; sizes.NumDiscStates  = 0;
sizes.NumOutputs     = (4+3)*5; sizes.NumInputs      = 0;
sizes.DirFeedthrough = 1; sizes.NumSampleTimes = 1;  
sys = simsizes(sizes);

x0  = [];
str = [];
ts  = [0 0];
simStateCompliance = 'UnknownSimState'; 
        
  case 3
    %Bodies:
    %1 - SLP
    %2 - SIP
    %3 - E-Field+
    %4 - E-Field-
        
    FR=Bodies.System.FR;
    scale=Bodies.System.scale; %To allow visualization of small perturbations
    
    %HUB
    Theta=Ang_ord(ceil((t+1/FR)*FR),:,1); %Obtain the correct instant of time in the ordered vector
    mat_ref=angle2dcm(0,-pi/2,-pi/2);     %To go to the defined ref sys in the 3d animation
    mat=mat_ref'*angle2dcm(Theta(3),Theta(2),Theta(1))'*mat_ref;  %Rotation matrix
    B1_r=vrrotmat2vec(mat);
    B1_p=vrcoordm2vr([0;0;0]');
    B1=[B1_r,B1_p];
    
    
    %Body 2
    Theta=Ang_ord(ceil((t+1/FR)*FR),:,2); %Obtain the correct instant of time in the ordered vector
    mat_ref=angle2dcm(-pi/2*0,-pi/2,-pi/2*0);
    mat=angle2dcm(pi/2,0,0)*mat_ref'*angle2dcm(Theta(3),Theta(2),Theta(1))'*mat_ref;
    B2_r=vrrotmat2vec(mat);
    

    %Body 3
    Theta=Ang_ord(ceil((t+1/FR)*FR),:,3); %Obtain the correct instant of time in the ordered vector
    mat_ref=angle2dcm(-pi,-pi/2,0);
    mat2=angle2dcm(0,pi/2,0);
    mat=mat2'*angle2dcm(0,pi/2,pi/2)*mat_ref'*angle2dcm(Theta(3),Theta(2),Theta(1))'*mat_ref*mat2;
    B3_r=vrrotmat2vec(mat);
    
    %Body 4
    Theta=Ang_ord(ceil((t+1/FR)*FR),:,4); %Obtain the correct instant of time in the ordered vector
    mat_ref=angle2dcm(-pi/2*0,-pi/2,-pi/2*0);
    mat=angle2dcm(pi/2,0,0)*mat_ref'*angle2dcm(Theta(3),Theta(2),Theta(1))'*mat_ref;
    B4_r=vrrotmat2vec(mat);
    
    %Body 5
    Theta=Ang_ord(ceil((t+1/FR)*FR),:,5); %Obtain the correct instant of time in the ordered vector
    mat_ref=angle2dcm(-pi/2*0,-pi/2,-pi/2*0);
    mat=angle2dcm(pi/2,0,0)*mat_ref'*angle2dcm(Theta(3),Theta(2),Theta(1))'*mat_ref;
    B5_r=vrrotmat2vec(mat);
    
    
    
    %Bodies translation
    B2_p=PHUB_3D_ord(ceil((t+1/FR)*FR),:,2)'+[0.11317;-0.16764;0.05448];
    B2_p=angle2dcm(pi/2,0,0)*B2_p+[0;0;0.15];
    B2=[B2_r,vrcoordm2vr(B2_p')];
    B3_p=PHUB_3D_ord(ceil((t+1/FR)*FR),:,3)'+[0.11317;-0.16764;0.05448];
    B3_p=angle2dcm(pi/2,0,0)*B3_p+[0;0;-0.15];
    B3=[B3_r,vrcoordm2vr(B3_p')];
    B4_p=PHUB_3D_ord(ceil((t+1/FR)*FR),:,4)'+[0.11317;-0.16764;0.05448];
    B4_p=angle2dcm(pi/2,0,0)*B4_p+[0;0;0.15];
    B4=[B4_r,vrcoordm2vr(B4_p')];   
    B5_p=PHUB_3D_ord(ceil((t+1/FR)*FR),:,5)'+[0.11317;-0.16764;0.05448];
    B5_p=angle2dcm(pi/2,0,0)*B5_p+[0;0;0.15];
    B5=[B5_r,vrcoordm2vr(B5_p')];
    
    
    sys=[B1,B2,B3,B4,B5];
end