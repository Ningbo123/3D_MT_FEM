 function Bdirc=get_bnd(Nx,Ny,Nz,NP)
  Bdirc=[];
   idxtop=1:(Ny+1)*(Nx+1);
   Bdirc=[Bdirc idxtop];
   
   idytop=idxtop+NP;
   Bdirc=[Bdirc idytop];
   
   idztop=idxtop+2*NP; 
   Bdirc=[Bdirc idztop];  
   
   iphitop=idxtop+3*NP;
   Bdirc=[Bdirc iphitop];
   
  for i=1:Nz-1
      idxmiddle=[i*(Ny+1)*(Nx+1)+(1:(Nx+1)) i*(Ny+1)*(Nx+1)+Ny*(Nx+1)+(1:(Nx+1))];     
      Bdirc=[Bdirc idxmiddle];
  
      
      idymiddle=idxmiddle+NP;
      Bdirc=[Bdirc idymiddle];
      idzmiddle=idxmiddle+2*NP; 
      Bdirc=[Bdirc idzmiddle];   
   
      iphimiddle=idxmiddle+3*NP;
      Bdirc=[Bdirc iphimiddle];       
  end
   
   for i=1:Nz-1
      for j=2:Ny
          idxmiddle=[i*(Ny+1)*(Nx+1)+(j-1)*(Nx+1)+1 i*(Ny+1)*(Nx+1)+(j-1)*(Nx+1)+Nx+1];
          Bdirc=[Bdirc idxmiddle];
          
          idymiddle=idxmiddle+NP;
          Bdirc=[Bdirc idymiddle];
   
          idzmiddle=idxmiddle+2*NP;
          Bdirc=[Bdirc idzmiddle]; 
   
         iphimiddle=idxmiddle+3*NP;
         Bdirc=[Bdirc iphimiddle];     
      end
   end  
   
   
  idxbot=Nz*(Ny+1)*(Nx+1)+1:(Nx+1)*(Ny+1)*(Nz+1);
  Bdirc=[Bdirc idxbot];
 
  idybot=idxbot+NP;
  Bdirc=[Bdirc idybot];
   
  idzbot=idxbot+2*NP; 
  Bdirc=[Bdirc idzbot];  
   
  iphibot=idxbot+3*NP;
  Bdirc=[Bdirc iphibot];
  Bdirc=Bdirc';
  Bdirc(:,2)=0;
  
   
   
   
   
   