function [ResXX,ResXY,ResYX,ResYY]=get_data(mu,f,NX,NY,Exx1,Eyy1,Hxx1,Hyy1,Exx2,Eyy2,Hxx2,Hyy2)
Zxx=zeros(NX,NY);
Zxy=zeros(NX,NY);
Zyx=zeros(NX,NY);
Zyy=zeros(NX,NY);
ResXX=zeros(NX,NY);
ResXY=zeros(NX,NY);
ResYX=zeros(NX,NY);
ResYY=zeros(NX,NY);
for i=1:NX
    for j=1:NY
        Temp=Hxx1(i,j)*Hyy2(i,j)-Hyy1(i,j)*Hxx2(i,j);
        Zxx(i,j)=( Exx1(i,j)*Hyy2(i,j)- Exx2(i,j)*Hyy1(i,j))/Temp;
        Zxy(i,j)=( Exx2(i,j)*Hxx1(i,j)- Exx1(i,j)*Hxx2(i,j))/Temp;
        Zyx(i,j)=( Eyy1(i,j)*Hyy2(i,j)-Eyy2(i,j)*Hyy1(i,j))/Temp;
        Zyy(i,j)=(Eyy2(i,j)*Hxx1(i,j)-Eyy1(i,j)*Hxx2(i,j))/Temp;
        ResYX(i,j)=abs((Zyx(i,j))^2*sqrt(-1)/(2*pi*f*mu));
        ResXY(i,j)=abs((Zxy(i,j))^2*sqrt(-1)/(2*pi*f*mu));
        ResXX(i,j)=abs((Zxx(i,j))^2*sqrt(-1)/(2*pi*f*mu));
        ResYY(i,j)=abs((Zyy(i,j))^2*sqrt(-1)/(2*pi*f*mu));
    end
end