
function [Ex1,Ey1,Ex2,Ey2,Ep_Hp1 Ep_Hp2]=MT1Danisinv(f,Z,Nair)


NZ=size(Z,2);
NE=NZ;      
NP=1+NE*3; 
for i=1:10
    rho(i,1)=10^6;
    rho(i,2)=10^6;
    rho(i,3)=10^6;
    alpha_S(i)=0;
    alpha_D(i)=0;
    alpha_L(i)=0;
end

for i=11:30
    rho(i,1)=100;
    rho(i,2)=100;
    rho(i,3)=100;
    alpha_S(i)=0;
    alpha_D(i)=0;
    alpha_L(i)=0;
end
for i=31:42
    rho(i,1)=50;
    rho(i,2)=50;
    rho(i,3)=50;
    alpha_S(i)=0;
    alpha_D(i)=0;
    alpha_L(i)=0;
end
for i=43:58
    rho(i,1)=200;
    rho(i,2)=200;
    rho(i,3)=200;
    alpha_S(i)=0;
    alpha_D(i)=0;
    alpha_L(i)=0;
    
end


for i=1:NE
    ME(1,i)=1+(i-1)*3;
    ME(2,i)=2+(i-1)*3;
    ME(3,i)=3+(i-1)*3;
    ME(4,i)=4+(i-1)*3;
end
for ff=1:size(f,2)
    K1=sparse(2*NP,2*NP);K2=sparse(2*NP,2*NP);K3=sparse(2*NP,2*NP);P=sparse(2*NP,1);
    for h=1:NE
        l=Z(h); 
        K=[148 -189   54 -13;-189  432 -297  54;54 -297  432 -189;-13 54 -189 148];
        for j=1:4
            NJ=ME(j,h);
            for k=1:4
                NK=ME(k,h);
                K1(NJ,NK)=K1(NJ,NK)+K(j,k)/40/l;
                K1(NP+NJ,NP+NK)=K1(NP+NJ,NP+NK)+K(j,k)/40/l;
            end
        end
    end
    %%%%%%%%%%%K2e%%%%%%%%%
    mu=(4e-7)*pi;
    w=2*pi*f(ff);
    m=sqrt(-1)*w*mu;
    for h=1:NE
        l=Z(h);
        K=[128 99 -36 19; 99 648 -81 -36; -36 -81 648 99; 19 -36 99 128];
        % ******************** ************************  %
        sigma1=1/rho(h,1);  
        sigma2=1/rho(h,2);  
        sigma3=1/rho(h,3);  
        sigma_zhuzhou=[sigma1 0 0; 0 sigma2 0; 0 0 sigma3]; 
        RzfuAlpha_S=[cos(-alpha_S(h)) sin(-alpha_S(h)) 0; -sin(-alpha_S(h)) cos(-alpha_S(h)) 0; 0 0 1]; 
        RxfuAlpha_D=[1 0 0; 0 cos(-alpha_D(h)) sin(-alpha_D(h)); 0 -sin(-alpha_D(h)) cos(-alpha_D(h))]; 
        RzfuAlpha_L=[cos(-alpha_L(h)) sin(-alpha_L(h)) 0; -sin(-alpha_L(h)) cos(-alpha_L(h)) 0; 0 0 1]; 
        RzAlpha_L=[cos(alpha_L(h)) sin(alpha_L(h)) 0; -sin(alpha_L(h)) cos(alpha_L(h)) 0; 0 0 1]; 
        RxAlpha_D=[1 0 0; 0 cos(alpha_D(h)) sin(alpha_D(h)); 0 -sin(alpha_D(h)) cos(alpha_D(h))]; 
        RzAlpha_S=[cos(alpha_S(h)) sin(alpha_S(h)) 0; -sin(alpha_S(h)) cos(alpha_S(h)) 0; 0 0 1]; 
        sigma_tensor=RzfuAlpha_S*RxfuAlpha_D*RzfuAlpha_L*sigma_zhuzhou*RzAlpha_L*RxAlpha_D*RzAlpha_S; 
        sigmaxx=sigma_tensor(1,1); 
        sigmaxy=sigma_tensor(1,2); 
        sigmaxz=sigma_tensor(1,3); 
        sigmayx=sigma_tensor(2,1); 
        sigmayy=sigma_tensor(2,2); 
        sigmayz=sigma_tensor(2,3); 
        sigmazx=sigma_tensor(3,1); 
        sigmazy=sigma_tensor(3,2); 
        sigmazz=sigma_tensor(3,3); 
        Axx=sigmaxx-sigmaxz*sigmazx/sigmazz;  
        Axy=sigmaxy-sigmaxz*sigmazy/sigmazz;  
        Ayx=sigmayx-sigmayz*sigmazx/sigmazz;  
        Ayy=sigmayy-sigmayz*sigmazy/sigmazz;  
        for j=1:4
            NJ=ME(j,h);
            for k=1:4
                NK=ME(k,h);
                K2(NJ,NK)=K2(NJ,NK)+K(j,k)*m*l*(Axx)/1680;
                K2(NJ,NP+NK)=K2(NJ,NP+NK)+K(j,k)*m*l*(Axy)/1680;
                K2(NP+NJ,NK)=K2(NP+NJ,NK)+K(j,k)*m*l*(Ayx)/1680;
                K2(NP+NJ,NP+NK)=K2(NP+NJ,NP+NK)+K(j,k)*m*l*(Ayy)/1680;
            end
        end
    end
    %%%%%%%%%%%%%%%%%K3e%%%%%%%%%%%%%%%%%%
    mu=(4e-7)*pi;
    w=2*pi*f(ff);
    m=sqrt(-1)*w*mu;
    a=sqrt(-m/2*(Axx+Ayy-sqrt((Axx-Ayy)^2+4*Axy*Ayx)));
    K3(NP,NP)=a;
    a=sqrt(-m/2*(Axx+Ayy+sqrt((Axx-Ayy)^2+4*Axy*Ayx)));
    K3(2*NP,2*NP)=a;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    v=K1-K2+K3;
    %********************** source 1  ******************************%
    DJ1=1; %        Ex
    v1=v;
    P1=P;
    v1(DJ1,DJ1)=v1(DJ1,DJ1)*10^10;
    P1(DJ1,1)=v1(DJ1,DJ1)*1;   %
    DJ2=NP+1;%        Ey
    v1(DJ2,DJ2)=v1(DJ2,DJ2)*10^10;
    P1(DJ2,1)=v1(DJ2,DJ2)*0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    DJ1=3*Nair+1;
    DJ2=NP+3*Nair+1;
    E(:,ff)=v1\P1;   %
    Ex1=[E(ME(1,:));E(NP)];
    Ey1=[E(ME(1,:)+NP);E(2*NP)];
    DExz(ff)=(-11*E(DJ1,ff)+18*E(DJ1+1,ff)-9*E(DJ1+2,ff)+2*E(DJ1+3,ff))/Z(11)/2;
    DEyz(ff)=(-11*E(DJ2,ff)+18*E(DJ2+1,ff)-9*E(DJ2+2,ff)+2*E(DJ2+3,ff))/Z(11)/2;
    ExA(ff)=E(DJ1,ff);
    HyA(ff)=DExz(ff)/(sqrt(-1)*w*mu);
    EyA(ff)=E(DJ2,ff);
    HxA(ff)=-DEyz(ff)/(sqrt(-1)*w*mu);
    Ep_Hp1=[ExA EyA HxA HyA];
    %**************************** source 2  *********************************%
    DJ1=1; %         Ex
    v2=v;
    P2=P;
    v2(DJ1,DJ1)=v2(DJ1,DJ1)*10^10;
    P2(DJ1,1)=v2(DJ1,DJ1)*0;   %
    DJ2=NP+1;%?        Ey
    v2(DJ2,DJ2)=v2(DJ2,DJ2)*10^10;
    P2(DJ2,1)=v2(DJ2,DJ2)*1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    DJ1=3*Nair+1;
    DJ2=NP+3*Nair+1;
    E(:,ff)=v2\P2;   %E
    Ex2=[E(ME(1,:));E(NP)];
    Ey2=[E(ME(1,:)+NP);E(2*NP)];
    DExz(ff)=(-11*E(DJ1,ff)+18*E(DJ1+1,ff)-9*E(DJ1+2,ff)+2*E(DJ1+3,ff))/Z(11)/2;
    DEyz(ff)=(-11*E(DJ2,ff)+18*E(DJ2+1,ff)-9*E(DJ2+2,ff)+2*E(DJ2+3,ff))/Z(11)/2;
    ExB(ff)=E(DJ1,ff);
    HyB(ff)=DExz(ff)/(sqrt(-1)*w*mu);
    EyB(ff)=E(DJ2,ff);
    HxB(ff)=-DEyz(ff)/(sqrt(-1)*w*mu);
    Ep_Hp2=[ExB EyB HxB HyB];
end

