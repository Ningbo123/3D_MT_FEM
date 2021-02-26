function [Exx2,Eyy2,Hxx2,Hyy2]=get_x2(v2,vm,vn,EN2,w,mu,Nx,Ny,Nz,NP,Nair,NL,a_x,b_y,c_z,Bdirc,Ep_Hp2,XYZ)
b1=vm*EN2;
b2=vn*EN2;
b22=[b1;b2];   
p2=b22;% 
b22=b22-v2(:,Bdirc(:,1))*b22(Bdirc(:,1));
b22(Bdirc(:,1))=Bdirc(:,2);
v2(Bdirc(:,1),:)=0;
v2(:,Bdirc(:,1))=0;
ii=Bdirc(:,1);
Kbnd=sparse(ii,ii,ones(size(ii)),4*NP,4*NP);
v2=v2+Kbnd;
tsolve3=clock;
tol=1e-10;maxsteps=5000;
setup.type='nofill';

[L,U]=ilu(v2,setup);
[xs,flag,relres,iter2,resvec2]=qmr(v2,b22,tol,maxsteps,L,U);
res2=resvec2/norm(b22);
save res2 res2
tsolve4=clock;
 disp('solve X2 time')
 etime(tsolve4,tsolve3)
  for i=1:NP
        Ax2(XYZ(1,i),XYZ(2,i),XYZ(3,i))=xs(i);
        Ay2(XYZ(1,i),XYZ(2,i),XYZ(3,i))=xs(NP+i);
        Az2(XYZ(1,i),XYZ(2,i),XYZ(3,i))=xs(NP*2+i);
        PHI2(XYZ(1,i),XYZ(2,i),XYZ(3,i))=xs(NP*3+i);
  end
    for i=1:Nx %x     
        for j=1:Ny % y    
            Axz1(i,j)=(-11* Ax2(i,j,Nair+1)+18* Ax2(i,j,Nair+2)-9* Ax2(i,j,Nair+3)+2* Ax2(i,j,Nair+4))/(6*c_z(Nair+1));
            Axz2(i,j)=(-11* Ax2(i+1,j,Nair+1)+18*Ax2(i+1,j,Nair+2)-9* Ax2(i+1,j,Nair+3)+2* Ax2(i+1,j,Nair+4))/(6*c_z(Nair+1));
            Axz3(i,j)=(-11* Ax2(i,j+1,Nair+1)+18*Ax2(i,j+1,Nair+2)-9* Ax2(i,j+1,Nair+3)+2* Ax2(i,j+1,Nair+4))/(6*c_z(Nair+1));
            Axz4(i,j)=(-11* Ax2(i+1,j+1,Nair+1)+18*Ax2(i+1,j+1,Nair+2)-9* Ax2(i+1,j+1,Nair+3)+2* Ax2(i+1,j+1,Nair+4))/(6*c_z(Nair+1));
            Ayz1(i,j)=(-11*Ay2(i,j,Nair+1)+18*Ay2(i,j,Nair+2)-9*Ay2(i,j,Nair+3)+2*Ay2(i,j,Nair+4))/(6*c_z(Nair+1));
            Ayz2(i,j)=(-11*Ay2(i+1,j,Nair+1)+18*Ay2(i+1,j,Nair+2)-9*Ay2(i+1,j,Nair+3)+2*Ay2(i+1,j,Nair+4))/(6*c_z(Nair+1));
            Ayz3(i,j)=(-11*Ay2(i,j+1,Nair+1)+18*Ay2(i,j+1,Nair+2)-9*Ay2(i,j+1,Nair+3)+2*Ay2(i,j+1,Nair+4))/(6*c_z(Nair+1));
            Ayz4(i,j)=(-11*Ay2(i+1,j+1,Nair+1)+18*Ay2(i+1,j+1,Nair+2)-9*Ay2(i+1,j+1,Nair+3)+2*Ay2(i+1,j+1,Nair+4))/(6*c_z(Nair+1));
            
            Ayz(i,j)=( Ayz1(i,j)+Ayz2(i,j)+Ayz3(i,j)+Ayz4(i,j))/4;
            Axz(i,j)=( Axz1(i,j)+ Axz2(i,j)+ Axz3(i,j)+ Axz4(i,j))/4;
            Azx(i,j)=((Az2(i+1,j,Nair+1)-Az2(i,j,Nair+1))/a_x(i)+(Az2(i+1,j+1,Nair+1)-Az2(i,j+1,Nair+1))/a_x(i))/2;
           
            Azy(i,j)=((Az2(i,j+1,Nair+1)-Az2(i,j,Nair+1))/b_y(j)+(Az2(i+1,j+1,Nair+1)-Az2(i+1,j,Nair+1))/b_y(j))/2;
         
            Exx2s(i,j)=(sqrt(-1)*w).*(( Ax2(i,j,Nair+1)+ Ax2(i+1,j,Nair+1)+ Ax2(i,j+1,Nair+1)+ Ax2(i+1,j+1,Nair+1))/4+((PHI2(i+1,j,Nair+1)-PHI2(i,j,Nair+1))/a_x(i)+(PHI2(i+1,j+1,Nair+1)-PHI2(i,j+1,Nair+1))/a_x(i))/2); %     Ĵų ȡΪ  Ԫ ĸ  ڵ ų   ƽ  ֵ
            Eyy2s(i,j)=(sqrt(-1)*w).*(( Ay2(i,j,Nair+1)+ Ay2(i+1,j,Nair+1)+ Ay2(i,j+1,Nair+1)+ Ay2(i+1,j+1,Nair+1))/4+((PHI2(i,j+1,Nair+1)-PHI2(i,j,Nair+1))/b_y(j)+(PHI2(i+1,j+1,Nair+1)-PHI2(i+1,j,Nair+1))/b_y(j))/2); %     Ĵų ȡΪ  Ԫ ĸ  ڵ ų   ƽ  ֵ

            Hyy2s(i,j)=(Axz(i,j)-Azx(i,j))/((mu));%  
            Hxx2s(i,j)=(Azy(i,j)-Ayz(i,j))/((mu));
        end
    end
    
    Exx2=Exx2s+Ep_Hp2(1);
    Eyy2=Eyy2s+Ep_Hp2(2);
    Hxx2=Hxx2s+Ep_Hp2(3);
    Hyy2=Hyy2s+Ep_Hp2(4);

end