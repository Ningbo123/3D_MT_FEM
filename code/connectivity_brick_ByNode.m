function [t] = connectivity_brick_ByNode(X,Y,Z)
nx = size(X,1);
ny = size(X,2);
nz = size(X,3);
Ne = (nx-1)*(ny-1)*(nz-1);% number of element

xn = X(:);%the coordinates of nodes
yn = Y(:);
zn = Z(:); 


p = [xn,yn,zn];%the xyz coordinate of nodes
t = zeros(Ne,8);
Num = reshape(1:size(p,1),nx,ny,nz); %Similary as X,Y,Z, but store the corresponding node number for X,Y,Z


%% Get the first node for each bricks which corresponding min(x),min(y),min(z) in each brick
NumTmp = Num(1:end-1,1:end-1,1:end-1);
t(:,1) = NumTmp(:);

%% Get the second node for each bricks which corresponding max(x),min(y),min(z) in each brick
NumTmp = Num(2:end,1:end-1,1:end-1);
t(:,2) = NumTmp(:);


%% Get the third node for each bricks which corresponding max(x),max(y),min(z) in each brick
NumTmp = Num(2:end,2:end,1:end-1);
t(:,3) = NumTmp(:);


%% Get the fourth node for each bricks which corresponding min(x),max(y),min(z) in each brick
NumTmp = Num(1:end-1,2:end,1:end-1);
t(:,4) = NumTmp(:);


%% Get the fifth node for each bricks which corresponding min(x),min(y),max(z) in each brick
NumTmp = Num(1:end-1,1:end-1,2:end);
t(:,5) = NumTmp(:);

%% Get the sixth node for each bricks which corresponding max(x),min(y),max(z) in each brick
NumTmp = Num(2:end,1:end-1,2:end);
t(:,6) = NumTmp(:);


%% Get the seventh node for each bricks which corresponding max(x),max(y),max(z) in each brick
NumTmp = Num(2:end,2:end,2:end);
t(:,7) = NumTmp(:);


%% Get the eight node for each bricks which corresponding min(x),max(y),max(z) in each brick
NumTmp = Num(1:end-1,2:end,2:end);
t(:,8) = NumTmp(:);
t=t';


