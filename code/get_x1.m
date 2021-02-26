function [Exx1,Eyy1,Hxx1,Hyy1]=get_x1(v1,vm,vn,EN1,w,mu,Nx,Ny,Nz,NP,Nair,NL,a_x,b_y,c_z,Bdirc,Ep_Hp1,XYZ)
b1=vm*EN1;
b2=vn*EN1;
b11=[b1;b2];
    
p1=b11;% 
b11=b11-v1(:,Bdirc(:,1))*b11(Bdirc(:,1));
b11(Bdirc(:,1))=Bdirc(:,2);
v1(Bdirc(:,1),:)=0;
v1(:,Bdirc(:,1))=0;
ii=Bdirc(:,1);
Kbnd=sparse(ii,ii,ones(size(ii)),4*NP,4*NP);
v1=v1+Kbnd;
%%    ILU+QMR
tol=1e-10;maxsteps=5000;
tsolve1=clock;
setup.type='nofill';
[L,U]=ilu(v1,setup); 
[xs,flag,relres,iter1,resvec1]=qmr(v1,b11,tol,maxsteps,L,U);
res1=resvec1/norm(b11);
save res1 res1
tsolve2=clock;
disp('solve X1 time')
etime(tsolve2,tsolve1)

save res1 res1

clear v1
  for i=1:NP
        Ax1(XYZ(1,i),XYZ(2,i),XYZ(3,i))=xs(i);
        Ay1(XYZ(1,i),XYZ(2,i),XYZ(3,i))=xs(NP+i);
        Az1(XYZ(1,i),XYZ(2,i),XYZ(3,i))=xs(NP*2+i);
        PHI1(XYZ(1,i),XYZ(2,i),XYZ(3,i))=xs(NP*3+i);
  end
  
  
    for i=1:Nx %x     
        for j=1:Ny % y    
            Axz1(i,j)=(-11* Ax1(i,j,Nair+1)+18* Ax1(i,j,Nair+2)-9* Ax1(i,j,Nair+3)+2* Ax1(i,j,Nair+4))/(6*c_z(Nair+1));
            Axz2(i,j)=(-11* Ax1(i+1,j,Nair+1)+18* Ax1(i+1,j,Nair+2)-9* Ax1(i+1,j,Nair+3)+2* Ax1(i+1,j,Nair+4))/(6*c_z(Nair+1));
            Axz3(i,j)=(-11* Ax1(i,j+1,Nair+1)+18* Ax1(i,j+1,Nair+2)-9* Ax1(i,j+1,Nair+3)+2* Ax1(i,j+1,Nair+4))/(6*c_z(Nair+1));
            Axz4(i,j)=(-11* Ax1(i+1,j+1,Nair+1)+18* Ax1(i+1,j+1,Nair+2)-9* Ax1(i+1,j+1,Nair+3)+2* Ax1(i+1,j+1,Nair+4))/(6*c_z(Nair+1));
            Ayz1(i,j)=(-11*Ay1(i,j,Nair+1)+18*Ay1(i,j,Nair+2)-9*Ay1(i,j,Nair+3)+2*Ay1(i,j,Nair+4))/(6*c_z(Nair+1));
            Ayz2(i,j)=(-11*Ay1(i+1,j,Nair+1)+18*Ay1(i+1,j,Nair+2)-9*Ay1(i+1,j,Nair+3)+2*Ay1(i+1,j,Nair+4))/(6*c_z(Nair+1));
            Ayz3(i,j)=(-11*Ay1(i,j+1,Nair+1)+18*Ay1(i,j+1,Nair+2)-9*Ay1(i,j+1,Nair+3)+2*Ay1(i,j+1,Nair+4))/(6*c_z(Nair+1));
            Ayz4(i,j)=(-11*Ay1(i+1,j+1,Nair+1)+18*Ay1(i+1,j+1,Nair+2)-9*Ay1(i+1,j+1,Nair+3)+2*Ay1(i+1,j+1,Nair+4))/(6*c_z(Nair+1));
            
            Ayz(i,j)=( Ayz1(i,j)+Ayz2(i,j)+Ayz3(i,j)+Ayz4(i,j))/4;
            Axz(i,j)=( Axz1(i,j)+ Axz2(i,j)+ Axz3(i,j)+ Axz4(i,j))/4;
            Azx(i,j)=((Az1(i+1,j,Nair+1)-Az1(i,j,Nair+1))/a_x(i)+(Az1(i+1,j+1,Nair+1)-Az1(i,j+1,Nair+1))/a_x(i))/2;
           
            Azy(i,j)=((Az1(i,j+1,Nair+1)-Az1(i,j,Nair+1))/b_y(j)+(Az1(i+1,j+1,Nair+1)-Az1(i+1,j,Nair+1))/b_y(j))/2;
         
            Exx1s(i,j)=(sqrt(-1)*w).*(( Ax1(i,j,Nair+1)+ Ax1(i+1,j,Nair+1)+ Ax1(i,j+1,Nair+1)+ Ax1(i+1,j+1,Nair+1))/4+((PHI1(i+1,j,Nair+1)-PHI1(i,j,Nair+1))/a_x(i)+(PHI1(i+1,j+1,Nair+1)-PHI1(i,j+1,Nair+1))/a_x(i))/2); %     Ĵų ȡΪ  Ԫ ĸ  ڵ ų   ƽ  ֵ
            Eyy1s(i,j)=(sqrt(-1)*w).*(( Ay1(i,j,Nair+1)+ Ay1(i+1,j,Nair+1)+ Ay1(i,j+1,Nair+1)+ Ay1(i+1,j+1,Nair+1))/4+((PHI1(i,j+1,Nair+1)-PHI1(i,j,Nair+1))/b_y(j)+(PHI1(i+1,j+1,Nair+1)-PHI1(i+1,j,Nair+1))/b_y(j))/2); %     Ĵų ȡΪ  Ԫ ĸ  ڵ ų   ƽ  ֵ

            Hyy1s(i,j)=(Axz(i,j)-Azx(i,j))/((mu));%  
            Hxx1s(i,j)=(Azy(i,j)-Ayz(i,j))/((mu));

        end
    end
    Exx1=Exx1s+Ep_Hp1(1);
    Eyy1=Eyy1s+Ep_Hp1(2);
    Hxx1=Hxx1s+Ep_Hp1(3);
    Hyy1=Hyy1s+Ep_Hp1(4);


  