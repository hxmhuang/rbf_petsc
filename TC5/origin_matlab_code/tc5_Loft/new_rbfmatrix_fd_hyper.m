function [DPx,DPy,DPz,L] = new_rbfmatrix_fd_hyper(atm,ep,fdsize,order,dim)
%%% [D,L] = rbfmatrix_fd_tree(x,ep,alpha,fdsize,order,dim)
% Requires kd-tree code.
% IN:
% x - nodes in this subroutine, NOT coordinate 'x' 
% ep - shape parameter
% fdsize - stencil size (good choices: 31, 50, 74, 101)
% order - L = L^order
% dim - dimension of Laplacian formula
% OUT:
% DPx - sparse differentiation matrix
% DPy - sparse differentiation matrix
% DPy - sparse differentiation matrix
% L - sparse dissipation matrix
x= [atm.pts.x atm.pts.y atm.pts.z];
N = length(x);

rbf = @(ep,rd2) exp(-ep^2*rd2);
drbf = @(ep,rd2) -2*ep^2*exp(-ep^2*rd2);

idx = knnsearch(x,x,'k',fdsize);
%dlmwrite('nn.md002.00009.txt',idx-1,' ')
MaskMatrix=zeros(N,N);

for j=1:N
    imat = idx(j,:);
    MaskMatrix(j,imat)=1;
end

dp1 = (-atm.pts.p_u * x' .* MaskMatrix)';
dp2 = (-atm.pts.p_v * x' .* MaskMatrix)';
dp3 = (-atm.pts.p_w * x' .* MaskMatrix)';

rd2= max(0,2*(1-x*x.'));
rbf_rd2 = rbf(ep,rd2);

%rd2v = drbf(ep,rd2);

rd2v =  -2*ep^2*rbf_rd2;

for j=1:N
    K0=rd2v(:,j);
    RightV1 = dp1(:,j).*K0;
    RightV2 = dp2(:,j).*K0;
    RightV3 = dp3(:,j).*K0;

    K1= repmat(MaskMatrix(j,:)',1,N);
    K2= repmat(MaskMatrix(j,:),N,1);
    K3= K1.*K2;
    K4= K3|eye(N);
    K5= rbf_rd2.* K4;
    ExtM=[K5 MaskMatrix(j,:)'; MaskMatrix(j,:) 0]
    
    ExtV=[RightV1;0];
    weights= ExtM\ExtV;
    DPx(:,j)= weights;
    
    ExtV=[RightV2;0]
    weights= ExtM\ExtV;
    DPy(:,j)= weights;

    ExtV=[RightV3;0]
    weights= ExtM\ExtV;
    DPz(:,j)= weights;
 
    
    K6=rd2(:,j).*MaskMatrix(j,:)';
    RightV4=ep^(2*order)*hyper(ep^2*K6,dim,order).*exp(-ep^2*K6).*MaskMatrix(j,:)';
    ExtV=[RightV4;0]
    weights= ExtM\ExtV;
    L(:,j)  = weights;
end
DPx=DPx(1:N,:)';
DPy=DPy(1:N,:)';
DPz=DPz(1:N,:)';
L=L';



function p=hyper(ep2r2,d,k)
%%% laplacian to the power k of dimension d
% ep2r2 - ep^2*r^2
n = length(ep2r2);
P = zeros(n,k+1);
P(:,1) = 1; P(:,2) = 4*ep2r2-2*d;
for j=3:k+1
    P(:,j) = 4*(ep2r2-2*j-d/2+4).*P(:,j-1) - 8*(j-2)*(2*j+d-6)*P(:,j-2);
end
p = P(:,k+1);
