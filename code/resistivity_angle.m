function [rho,rho0,alpha_S,alpha_S0,alpha_D,alpha_D0,alpha_L,alpha_L0] =resistivity_angle(Nair,NX,NY,NZ)


h=1:Nair*NX*NY;
rho(h,1)=1000000.d0;rho(h,2)=1000000.d0;rho(h,3)=1000000.d0;
alpha_S(h)=0.d0;alpha_D(h)=0.d0;alpha_L(h)=0.d0;

alpha_S0(h)=0.d0;alpha_D0(h)=0.d0;alpha_L0(h)=0.d0;
rho0(h,1)=1000000.d0;rho0(h,2)=1000000.d0;rho0(h,3)=1000000.d0;%
for h=Nair*NX*NY+1:NZ*NX*NY
rho(h,1)=100.d0;rho(h,2)=100.d0;rho(h,3)=100.d0;
rho0(h,1)=100.d0;rho0(h,2)=100.d0;rho0(h,3)=100.d0;
alpha_S(h)=0.d0;alpha_D(h)=0.d0;alpha_L(h)=0.d0;
alpha_S0(h)=0.d0;alpha_D0(h)=0.d0;alpha_L0(h)=0.d0;
end
%% the resistivity and Eular's angles of the anomalies
for L=17:30 
    for M=29:42
        for N=25:28
            h=(L-1)*NX*NY+(M-1)*NX+N;
            rho(h,1)=500.d0;rho(h,2)=50.d0;rho(h,3)=30.d0;    
            alpha_S(h)=pi/6;alpha_D(h)=pi/4;alpha_L(h)=pi/3;
        end
    end
end
for L=17:30 
    for M=29:42
        for N=43:46
            h=(L-1)*NX*NY+(M-1)*NX+N;
            rho(h,1)=500.d0;rho(h,2)=50.d0;rho(h,3)=30.d0;  
            rho0(h,1)=100.d0;rho0(h,2)=100.d0;rho0(h,3)=100.d0;  
            alpha_S(h)=pi/6;alpha_D(h)=pi/4;alpha_L(h)=pi/3;
            alpha_S0(h)=0.d0;alpha_D0(h)=0.d0;alpha_L0(h)=0.d0;
        end
    end
end