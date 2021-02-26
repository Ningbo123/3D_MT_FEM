function [EN1,EN2]=get_Ex0(Ex1,Ex2,Ey1,Ey2,Nx,Ny,Nz)
DZ=size(Ex1,1);
EX1=Ex1;EX2=Ex2;
EY1=Ey1;EY2=Ey2;
EZ1=zeros(DZ,1);EZ2=zeros(DZ,1);
   DC=(Nx+1)*(Ny+1);
   index=reshape(repmat(1:Nz+1,DC,1),DC*(Nz+1),1);
   Exx1=EX1(index);
   Eyy1=EY1(index);
   Ezz1=EZ1(index);
   EN1=[Exx1;Eyy1;Ezz1].';
   Exx2=EX2(index);
   Eyy2=EY2(index);
   Ezz2=EZ2(index);
   EN2=[Exx2;Eyy2;Ezz2].';
   EN1=EN1.';
   EN2=EN2.';
   
   
   
   
   
   