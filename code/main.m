clear all
clc
tic
format long
Input;
global mu
mu=(4e-7)*pi;
epsilo0=8.8542*(10^(-12));
f=0.001;
w=2*pi*f;
NX=length(A_X);
NY=length(B_Y);
NZ=length(C_Z);
NE=NX*NY*NZ;%the number of cells
NP=(NX+1)*(NY+1)*(NZ+1);%the number of nodes in the domain
NL=(NX+1)*NY*NZ+NX*(NY+1)*NZ+NX*NY*(NZ+1);
i=0;
for ll=1:NZ+1
    for mm=1:NY+1
        for nn=1:NX+1
            i=i+1;
            XYZ(1:3,i)=[nn mm ll]  ;  
        end
    end
end

fprintf('%s \n','Stage 1: Prepare for nodes and connectivity matrix');
[X,Y,Z]=ndgrid(x,y,z);
ME= connectivity_brick_ByNode(X,Y,Z);
[rho,rho0,alpha_S,alpha_S0,alpha_D,alpha_D0,alpha_L,alpha_L0] =resistivity_angle(Nair,NX,NY,NZ);
fprintf('%s \n','Stage 2: Produce global stiffness matrix');
%% get the stiffnes matrix v and right hand corresponding matrix v1 v2
[v,v1,v2]=get_stiffness_matrix(A_X,B_Y,C_Z,NE,NX,NY,NZ,w,rho,rho0,alpha_S,alpha_D,alpha_L,alpha_S0,alpha_D0,alpha_L0,mu,ME,NP);%
clear rho rho0 alpha_S alpha_S0 alpha_D alpha_D0 alpha_L alpha_L0
fprintf('%s \n','Stage 3: add the boundary condition and solve the equation');
[Ex1,Ey1,Ex2,Ey2,Ep_Hp1,Ep_Hp2]=MT1Danisinv(f,C_Z,Nair);
Bdirc=get_bnd(NX,NY,NZ,NP);
[EN1,EN2]=get_Ex0(Ex1,Ex2,Ey1,Ey2,NX,NY,NZ);
[Exx1,Eyy1,Hxx1,Hyy1]=get_x1(v,v1,v2,EN1,w,mu,NX,NY,NZ,NP,Nair,NL,A_X,B_Y,C_Z,Bdirc,Ep_Hp1,XYZ);
[Exx2,Eyy2,Hxx2,Hyy2]=get_x2(v,v1,v2,EN2,w,mu,NX,NY,NZ,NP,Nair,NL,A_X,B_Y,C_Z,Bdirc,Ep_Hp2,XYZ);
fprintf('%s \n','Stage 4: get the apparent resistivity');
[ResXX,ResXY,ResYX,ResYY]=get_data(mu,f,NX,NY,Exx1,Eyy1,Hxx1,Hyy1,Exx2,Eyy2,Hxx2,Hyy2)
save ResXX ResXX
save ResXY ResXY
save ResYX ResYX
save ResYY ResYY
painting1;