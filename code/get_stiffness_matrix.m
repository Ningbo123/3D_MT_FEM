function [v,v1,v2]=get_stiffness_matrix(A_X,B_Y,C_Z,NE,NX,NY,NZ,w,rho,rho0,alpha_S,alpha_D,alpha_L,alpha_S0,alpha_D0,alpha_L0,mu,t,NP)

                   
k1=[4 -4 -2 2 2 -2 -1 1;
    -4 4 2 -2 -2  2 1 -1;
    -2 2 4 -4 -1  1 2 -2;
    2 -2 -4 4 1 -1 -2 2;
    2 -2 -1 1 4 -4 -2 2;
    -2 2 1 -1 -4 4 2 -2;
    -1 1 2 -2 -2 2 4 -4;
    1 -1 -2 2 2 -2 -4 4];
k2=[4 2 -2 -4 2 1 -1 -2;
    2 4 -4 -2 1 2 -2 -1;
    -2 -4 4 2 -1 -2 2 1;
    -4 -2 2 4 -2 -1 1 2;
    2 1 -1 -2 4 2 -2 -4;
    1 2 -2 -1 2 4 -4 -2;
    -1 -2 2 1 -2 -4 4 2;
    -2 -1 1 2 -4 -2 2 4];
k3=[4 2 1 2 -4 -2 -1 -2;
    2 4 2 1 -2 -4 -2 -1;
    1 2 4 2 -1 -2 -4 -2;
    2 1 2 4 -2 -1 -2 -4;
    -4 -2 -1 -2 4 2 1 2;
    -2 -4 -2 -1 2 4 2 1;
    -1 -2 -4 -2 1 2 4 2;
    -2 -1 -2 -4 2 1 2 4];
k4=[8 4 2 4 4 2 1 2;
    4 8 4 2 2 4 2 1;
    2 4 8 4 1 2 4 2;
    4 2 4 8 2 1 2 4;
    4 2 1 2 8 4 2 4;
    2 4 2 1 4 8 4 2;
    1 2 4 2 2 4 8 4;
    2 1 2 4 4 2 4 8];
%%
k5=[-4 4 2 -2 -2 2 1 -1;
    -4 4 2 -2 -2 2 1 -1;
    -2 2 4 -4 -1 1 2 -2;
    -2 2 4 -4 -1 1 2 -2;
    -2 2 1 -1 -4 4 2 -2;
    -2 2 1 -1 -4 4 2 -2;
    -1 1 2 -2 -2 2 4 -4;
    -1 1 2 -2 -2 2 4 -4];%
k6=[-4 -2 2 4 -2 -1 1 2;
    -2 -4 4 2 -1 -2 2 1;
    -2 -4 4 2 -1 -2 2 1;
    -4 -2 2 4 -2 -1 1 2;
    -2 -1 1 2 -4 -2 2 4;
    -1 -2 2 1 -2 -4 4 2;
    -1 -2 2 1 -2 -4 4 2;
    -2 -1 1 2 -4 -2 2 4];%
k7=[-4 -2 -1 -2 4 2 1 2;
    -2 -4 -2 -1 2 4 2 1;
    -1 -2 -4 -2 1 2 4 2;
    -2 -1 -2 -4 2 1 2 4;
    -4 -2 -1 -2 4 2 1 2;
    -2 -4 -2 -1 2 4 2 1;
    -1 -2 -4 -2 1 2 4 2;
    -2 -1 -2 -4 2 1 2 4];%
%%
k8=[4 -4 -2 2 2 -2 -1 1;
    -4 4 2 -2 -2 2 1 -1;
    -2 2 4 -4 -1 1 2 -2;
    2 -2 -4 4 1 -1 -2 2;
    2 -2 -1 1 4 -4 -2 2;
    -2 2 1 -1 -4 4 2 -2;
    -1 1 2 -2 -2 2 4 -4
    1 -1 -2 2 2 -2 -4 4];
k9=[6 -6 -6 6 3 -3 -3 3;
    6 -6 -6 6 3 -3 -3 3;
    -6 6 6 -6 -3 3 3 -3;
    -6 6 6 -6 -3 3 3 -3;
    3 -3 -3 3 6 -6 -6 6;
    3 -3 -3 3 6 -6 -6 6;
    -3 3 3 -3 -6 6 6 -6;
    -3 3 3 -3 -6 6 6 -6];
k10=[6 -6 -3 3 6 -6 -3 3;
    6 -6 -3 3 6 -6 -3 3;
    3 -3 -6 6 3 -3 -6 6;
    3 -3 -6 6 3 -3 -6 6;
    -6 6 3 -3 -6 6 3 -3;
    -6 6 3 -3 -6 6 3 -3;
    -3 3 6 -6 -3 3 6 -6;
    -3 3 6 -6 -3 3 6 -6];
k11=[6 6 -6 -6 3 3 -3 -3;
    -6 -6 6 6 -3 -3 3 3;
    -6 -6 6 6 -3 -3 3 3;
    6 6 -6 -6 3 3 -3 -3;
    3 3 -3 -3 6 6 -6 -6;
    -3 -3 3 3 -6 -6 6 6;
    -3 -3 3 3 -6 -6 6 6;
    3 3 -3 -3 6 6 -6 -6];

k12=[4 2 -2 -4 2 1 -1 -2;
    2 4 -4 -2 1 2 -2 -1;
    -2 -4 4 2 -1 -2 2 1;
    -4 -2 2 4  -2 -1 1 2;
    2 1 -1 -2 4 2 -2 -4;
    1 2 -2 -1 2 4 -4 -2;
    -1 -2 2 1 -2 -4 4 2;
    -2 -1 1 2 -4 -2 2 4];
k13=[6 3 -3 -6 6 3 -3 -6;
    3 6 -6 -3 3 6 -6 -3;
    3 6 -6 -3 3 6 -6 -3;
    6 3 -3 -6 6 3 -3 -6;
    -6 -3 3 6 -6 -3 3 6;
    -3 -6 6 3 -3 -6 6 3;
    -3 -6 6 3 -3 -6 6 3;
    -6 -3 3 6 -6 -3 3 6];%
k14=[6 6 3 3 -6 -6 -3 -3;
    -6 -6 -3 -3 6 6 3 3;
    -3 -3 -6 -6 3 3 6 6;
    3 3 6 6 -3 -3 -6 -6;
    6 6 3 3 -6 -6 -3 -3;
    -6 -6 -3 -3 6 6 3 3;
    -3 -3 -6 -6 3 3 6 6;
    3 3 6 6 -3 -3 -6 -6];%
k15=[6 3 3 6 -6 -3 -3 -6;
    3 6 6 3 -3 -6 -6 -3;
    -3 -6 -6 -3 3 6 6 3;
    -6 -3 -3 -6 6 3 3 6;
    6 3 3 6 -6 -3 -3 -6;
    3 6 6 3 -3 -6 -6 -3;
    -3 -6 -6 -3 3 6 6 3;
    -6 -3 -3 -6 6 3 3 6];%
k16=[4 2 1 2 -4 -2 -1 -2;
    2 4 2 1 -2 -4 -2 -1;
    1 2 4 2 -1 -2 -4 -2;
    2 1 2 4 -2 -1 -2 -4;
    -4 -2 -1 -2 4 2 1 2;
    -2 -4 -2 -1 2 4 2 1;
    -1 -2 -4 -2 1 2 4 2;
    -2 -1 -2 -4 2 1 2 4];%

for h=1:NE 
  %% the principal conductivity tensor of the model and background
  sigma_principal=[1/rho(h,1) 0 0;
      0 1/rho(h,2) 0;
      0 0 1/rho(h,3)];
  sigma_principal_background=[1/rho0(h,1)      0                0; ...
      0             1/rho0(h,2)         0; ...
      0                0         1/rho0(h,3)];
  %% Rotation transformation matrix R(-aS),R(-aD),R(-aL) of the model and background
  RzfuAlpha_S=[cos(-alpha_S(h)) sin(-alpha_S(h))   0;
      -sin(-alpha_S(h)) cos(-alpha_S(h))      0;
      0                   0                 1];
  RzfuAlpha_S0=[cos(-alpha_S0(h)) sin(-alpha_S0(h))  0;
      -sin(-alpha_S0(h))  cos(-alpha_S0(h)) 0;
      0                     0             1];
  RxfuAlpha_D=[1          0                       0;
      0 cos(-alpha_D(h))  sin(-alpha_D(h));
      0 -sin(-alpha_D(h)) cos(-alpha_D(h))];
  RxfuAlpha_D0=[1 0 0;
      0 cos(-alpha_D0(h)) sin(-alpha_D0(h));
      0 -sin(-alpha_D0(h)) cos(-alpha_D0(h))];
  RzfuAlpha_L=[cos(-alpha_L(h)) sin(-alpha_L(h)) 0;
      -sin(-alpha_L(h)) cos(-alpha_L(h)) 0;
      0 0 1];
  RzfuAlpha_L0=[cos(-alpha_L0(h)) sin(-alpha_L0(h)) 0;
      -sin(-alpha_L0(h)) cos(-alpha_L0(h)) 0;
      0 0 1];
  %%  Rotation transformation matrix R(aS),R(aD),R-aL) of the model and background
  RzAlpha_L=[cos(alpha_L(h)) sin(alpha_L(h)) 0;
      -sin(alpha_L(h)) cos(alpha_L(h)) 0;
      0 0 1];
  RzAlpha_L0=[cos(alpha_L0(h)) sin(alpha_L0(h)) 0;
      -sin(alpha_L0(h)) cos(alpha_L0(h)) 0;
      0 0 1];
  
  RxAlpha_D=[1 0 0;
      0 cos(alpha_D(h)) sin(alpha_D(h));
      0 -sin(alpha_D(h)) cos(alpha_D(h))];
  RxAlpha_D0=[1 0 0;
      0 cos(alpha_D0(h)) sin(alpha_D0(h));
      0 -sin(alpha_D0(h)) cos(alpha_D0(h))];
  RzAlpha_S=[cos(alpha_S(h)) sin(alpha_S(h)) 0;
      -sin(alpha_S(h)) cos(alpha_S(h)) 0;
      0 0 1];
  RzAlpha_S0=[cos(alpha_S0(h)) sin(alpha_S0(h)) 0;
      -sin(alpha_S0(h)) cos(alpha_S0(h)) 0;
      0 0 1];
  %% the calculation of conductivity tensor sigma=R(-aS)*R(-aD)*R(-aL)*the principal conductivity tensor*R(aL)*R(aD)*R(aS)
  sigma_tensor=RzfuAlpha_S*RxfuAlpha_D*RzfuAlpha_L*sigma_principal*RzAlpha_L*RxAlpha_D*RzAlpha_S;
  sigma_tensorp=RzfuAlpha_S0*RxfuAlpha_D0*RzfuAlpha_L0*sigma_principal_background*RzAlpha_L0*RxAlpha_D0*RzAlpha_S0;%�絼�������ļ���
  %%
  sigma_xx=sigma_tensor(1,1);
  sigma_xy=sigma_tensor(1,2);
  sigma_xz=sigma_tensor(1,3);
  sigma_yx=sigma_tensor(2,1);
  sigma_yy=sigma_tensor(2,2);
  sigma_yz=sigma_tensor(2,3);
  sigma_zx=sigma_tensor(3,1);
  sigma_zy=sigma_tensor(3,2);
  sigma_zz=sigma_tensor(3,3);
  %% the calculation of sigma-sigmap (the model conductivity tensor-the background conductivity tensor)
  delsigma_xx=sigma_tensor(1,1)-sigma_tensorp(1,1);
  delsigma_xy=sigma_tensor(1,2)-sigma_tensorp(1,2);
  delsigma_xz=sigma_tensor(1,3)-sigma_tensorp(1,3);
  delsigma_yx=sigma_tensor(2,1)-sigma_tensorp(2,1);
  delsigma_yy=sigma_tensor(2,2)-sigma_tensorp(2,2);
  delsigma_yz=sigma_tensor(2,3)-sigma_tensorp(2,3);
  delsigma_zx=sigma_tensor(3,1)-sigma_tensorp(3,1);
  delsigma_zy=sigma_tensor(3,2)-sigma_tensorp(3,2);
  delsigma_zz=sigma_tensor(3,3)-sigma_tensorp(3,3);
  sigma_inv=inv(sigma_tensor);
%   a=mod(h-1,NX)+1;
%   b=mod(h-a_xn,NX*NY)/NX+1;
%   c=floor((h-1)/(NX*NY))+1;
a_xn=mod(h-1,NX)+1;
b_yn=mod(h-a_xn,NX*NY)/NX+1;
c_zn=floor((h-1)/(NX*NY))+1;
a=A_X(a_xn);
b=B_Y(b_yn);
c=C_Z(c_zn);
  s11factor=-b*c/36/a*k1-a*c/36/b*k2-a*b/36/c*k3+sqrt(-1)*w*mu*sigma_xx*a*b*c/216*k4;
  s12factor= sqrt(-1)*w*mu*sigma_xy*a*b*c/216*k4;
  s13factor= sqrt(-1)*w*mu*sigma_xz*a*b*c/216*k4;%
  s14factor= sqrt(-1)*w*mu*(b*c*sigma_xx/72*k5+a*c*sigma_xy/72*k6+a*b/72*sigma_xz*k7);%
  s21factor=sqrt(-1)*w*mu*sigma_yx*a*b*c/216*k4;%
  s22factor=-b*c/36/a*k1-a*c/36/b*k2-a*b/36/c*k3+sqrt(-1)*w*mu*sigma_yy*a*b*c/216*k4;
  s23factor=sqrt(-1)*w*mu*sigma_yz*a*b*c/216*k4;%
  s24factor=sqrt(-1)*w*mu*(b*c*sigma_yx/72*k5+a*c*sigma_yy/72*k6+a*b/72*sigma_yz*k7);%
  s31factor=sqrt(-1)*w*mu*sigma_zx*a*b*c/216*k4;%
  s32factor=sqrt(-1)*w*mu*sigma_zy*a*b*c/216*k4;%
  s33factor=-b*c/36/a*k1-a*c/36/b*k2-a*b/36/c*k3+sqrt(-1)*w*mu*sigma_zz*a*b*c/216*k4;
  s34factor=sqrt(-1)*w*mu*(b*c*sigma_zx/72*k5+a*c*sigma_zy/72*k6+a*b/72*sigma_zz*k7);%
  s41factor=s14factor.';
  s42factor=s24factor.';
  s43factor=s34factor.';
  s44factor=sqrt(-1)*w*mu*[b*c*sigma_xx/36/a*k8+c*sigma_yx/72*k9+sigma_zx*b/72*k10+ ...
      sigma_xy*c/72*k11+a*c*sigma_yy/36/b*k12+a*sigma_zy/72*k13+...
      +b*sigma_xz/72*k14+a*sigma_yz/72*k15+a*b*sigma_zz/36/c*k16];%
  
  y11factor=-mu*(delsigma_xx)*a*b*c/216*k4;
  y12factor=-mu*(delsigma_xy)*a*b*c/216*k4;
  y13factor=-mu*(delsigma_xz)*a*b*c/216*k4;
  y21factor=-mu*(delsigma_yx)*a*b*c/216*k4;
  y22factor=-mu*(delsigma_yy)*a*b*c/216*k4;
  y23factor=-mu*(delsigma_yz)*a*b*c/216*k4;
  y31factor=-mu*(delsigma_zx)*a*b*c/216*k4;
  y32factor=-mu*(delsigma_zy)*a*b*c/216*k4;
  y33factor=-mu*(delsigma_zz)*a*b*c/216*k4;
  y41factor=-mu*(b*c*delsigma_xx/72*k5.'+a*c*delsigma_yx/72*k6.'+a*b/72*delsigma_zx*k7.');
  y42factor=-mu*(b*c*delsigma_xy/72*k5.'+a*c*delsigma_yy/72*k6.'+a*b/72*delsigma_zy*k7.');
  y43factor=-mu*(b*c*delsigma_xz/72*k5.'+a*c*delsigma_yz/72*k6.'+a*b/72*delsigma_zz*k7.');
  s11e=reshape(s11factor.',8*8,1);
  s12e=reshape(s12factor.',8*8,1);
  s13e=reshape(s13factor.',8*8,1);
  s14e=reshape(s14factor.',8*8,1);
  s21e=reshape(s21factor.',8*8,1);
  s22e=reshape(s22factor.',8*8,1);
  s23e=reshape(s23factor.',8*8,1);
  s24e=reshape(s24factor.',8*8,1);
  s31e=reshape(s31factor.',8*8,1);
  s32e=reshape(s32factor.',8*8,1);
  s33e=reshape(s33factor.',8*8,1);
  s34e=reshape(s34factor.',8*8,1);
  
  s41e=reshape(s41factor.',8*8,1);
  s42e=reshape(s42factor.',8*8,1);
  s43e=reshape(s43factor.',8*8,1);
  s44e=reshape(s44factor.',8*8,1);
  y11=reshape(y11factor.',8*8,1);
  y12=reshape(y12factor.',8*8,1);
  y13=reshape(y13factor.',8*8,1);
  y21=reshape(y21factor.',8*8,1);
  y22=reshape(y22factor.',8*8,1);
  y23=reshape(y23factor.',8*8,1);
  y31=reshape(y31factor.',8*8,1);
  y32=reshape(y32factor.',8*8,1);
  y33=reshape(y33factor.',8*8,1);
  y4a=reshape(y41factor.',8*8,1);
  y4b=reshape(y42factor.',8*8,1);
  y4c=reshape(y43factor.',8*8,1);
  S11(h:NE:(8*8-1)*NE+h)=s11e;
  S12(h:NE:(8*8-1)*NE+h)=s12e;
  S13(h:NE:(8*8-1)*NE+h)=s13e;
  S14(h:NE:(8*8-1)*NE+h)=s14e;
  S21(h:NE:(8*8-1)*NE+h)=s21e;
  S22(h:NE:(8*8-1)*NE+h)=s22e;
  S23(h:NE:(8*8-1)*NE+h)=s23e;
  S24(h:NE:(8*8-1)*NE+h)=s24e;
  S31(h:NE:(8*8-1)*NE+h)=s31e;
  S32(h:NE:(8*8-1)*NE+h)=s32e;
  S33(h:NE:(8*8-1)*NE+h)=s33e;
  S34(h:NE:(8*8-1)*NE+h)=s34e;
  S41(h:NE:(8*8-1)*NE+h)=s41e;
  S42(h:NE:(8*8-1)*NE+h)=s42e;
  S43(h:NE:(8*8-1)*NE+h)=s43e;
  S44(h:NE:(8*8-1)*NE+h)=s44e;
  Y11(h:NE:(8*8-1)*NE+h)=y11;
  Y12(h:NE:(8*8-1)*NE+h)=y12;
  Y13(h:NE:(8*8-1)*NE+h)=y13;
  Y21(h:NE:(8*8-1)*NE+h)=y21;
  Y22(h:NE:(8*8-1)*NE+h)=y22;
  Y23(h:NE:(8*8-1)*NE+h)=y23;
  Y31(h:NE:(8*8-1)*NE+h)=y31;
  Y32(h:NE:(8*8-1)*NE+h)=y32;
  Y33(h:NE:(8*8-1)*NE+h)=y33;
  Y4a(h:NE:(8*8-1)*NE+h)=y4a;
  Y4b(h:NE:(8*8-1)*NE+h)=y4b;
  Y4c(h:NE:(8*8-1)*NE+h)=y4c;
end
ii=zeros(8*8*NE,1);
jj=zeros(8*8*NE,1);
index=0;
for it1=1:8
    for jt1=1:8
        ii(index+1:index+NE)=t(it1,:);
        jj(index+1:index+NE)=t(jt1,:);
        index=index+NE;
    end
end
v11=sparse(ii,jj,S11,NP,NP);
v12=sparse(ii,jj,S12,NP,NP);
v13=sparse(ii,jj,S13,NP,NP);
v14=sparse(ii,jj,S14,NP,NP);
v22=sparse(ii,jj,S22,NP,NP);
v23=sparse(ii,jj,S23,NP,NP);
v24=sparse(ii,jj,S24,NP,NP);
v33=sparse(ii,jj,S33,NP,NP);
v34=sparse(ii,jj,S34,NP,NP);
v44=sparse(ii,jj,S44,NP,NP);
clear S11 S12 S13 S14  S22 S23 S24 S33 S34 S44
v=[v11 v12 v13 v14; v12.' v22 v23 v24; v13.' v23.' v33 v34; v14.' v24.' v34.' v44];
clear v11 v12 v13 v14 v22 v23 v24 v33 v34 v44
Y11a=sparse(ii,jj,Y11,NP,NP);
Y12a=sparse(ii,jj,Y12,NP,NP);
Y13a=sparse(ii,jj,Y13,NP,NP);
Y21a=sparse(ii,jj,Y21,NP,NP);
Y22a=sparse(ii,jj,Y22,NP,NP);
Y23a=sparse(ii,jj,Y23,NP,NP);
Y31a=sparse(ii,jj,Y31,NP,NP);
Y32a=sparse(ii,jj,Y32,NP,NP);
Y33a=sparse(ii,jj,Y33,NP,NP);
v1=[Y11a Y12a Y13a;Y21a Y22a Y23a; Y31a Y32a Y33a];
clear Y11a Y12a Y13a Y21a Y22a Y23a Y31a Y32a Y33a
v2=[sparse(ii,jj,Y4a,NP,NP) sparse(ii,jj,Y4b,NP,NP) sparse(ii,jj,Y4c,NP,NP)];
clear Y4a Y4b Y4c

