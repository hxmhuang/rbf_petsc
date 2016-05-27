function [DPx,DPy,DPz,L] = rbfmatrix_fd_hyper(x,ep,fdsize,order,dim,Px,Py,Pz)
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

N = length(x);
srange = sqrt(6*fdsize/N); % search range

rbf = @(ep,rd2) exp(-ep^2*rd2);
drbf = @(ep,rd2) -2*ep^2*exp(-ep^2*rd2);

weightsDx = zeros(N*fdsize,1);
weightsDy = zeros(N*fdsize,1);
weightsDz = zeros(N*fdsize,1);
weightsL = zeros(N*fdsize,1);

ind_i = zeros(N*fdsize,1);
ind_j = zeros(N*fdsize,1);

A = ones(fdsize+1,fdsize+1); A(end,end) = 0;
B = zeros(fdsize+1,1);

idx = knnsearch(x,x,'k',fdsize);
MaskMatrix=zeros(N,fdsize);
A_new = ones(N+1,N+1); A(end,end) = 0;
B_new1 = zeros(N+1,1);
B_new2 = zeros(N+1,1);
B_new3 = zeros(N+1,1);

for j=1:N
    imat = idx(j,:);
    ind_i((j-1)*fdsize+1:j*fdsize) = j;
    ind_j((j-1)*fdsize+1:j*fdsize) = imat;
    MaskMatrix(j,imat)=1
end

%MaskMatrix = sparse(ind_i,ind_j,1,N,N); 

dp1 = -Px*x'.*MaskMatrix;
dp2 = -Py*x'.*MaskMatrix;
dp3 = -Pz*x'.*MaskMatrix;

dp1=dp1';
dp2=dp2';
dp3=dp3';

%rd2_new= max(0,2*(1-x*x').*MaskMatrix);
rd2_new= max(0,2*(1-x*x.'));
A_new = rbf(ep,rd2_new);

%rd2v_new = drbf(ep,rd2_new);

rd2v_new =  -2*ep^2*A_new;

for j=1:N
    M_tmp1=rd2v_new(:,j);
    RightV1 = dp1(:,j).*M_tmp1;
    RightV2 = dp2(:,j).*M_tmp1;
    RightV3 = dp3(:,j).*M_tmp1;

    K_tmp1= repmat(MaskMatrix(j,:)',1,N);
    K_tmp2= repmat(MaskMatrix(j,:),N,1);
    K_tmp3= K_tmp1.*K_tmp2;
    K_tmp4= K_tmp3|eye(N);
    K_tmp5= A_new.* K_tmp4;
    ExtM=[K_tmp5 MaskMatrix(j,:)'; MaskMatrix(j,:) 0]
    
    ExtV=[RightV1;0];
    weights1_new = ExtM\ExtV;
    DPx_new(j,:)= weights1_new(1:N);
    
    ExtV=[RightV2;0]
    weights2_new = ExtM\ExtV;
    DPy_new(j,:)= weights2_new(1:N);

    ExtV=[RightV3;0]
    weights3_new = ExtM\ExtV;
    DPz_new(j,:)= weights3_new(1:N);
 
    
    M_tmp2=rd2_new(:,j).*MaskMatrix(j,:)';
    RightV4=ep^(2*order)*hyper(ep^2*M_tmp2,dim,order).*exp(-ep^2*M_tmp2).*MaskMatrix(j,:)';
    ExtV=[RightV4;0]
    weights4_new = ExtM\ExtV;
    L_new(j,:)= weights4_new(1:N);
end


for j=1:N
    
    imat = idx(j,:);
    ind_i((j-1)*fdsize+1:j*fdsize) = j;
    ind_j((j-1)*fdsize+1:j*fdsize) = imat;
    
    dp_tmp = (x(j,1)*x(imat,1) + x(j,2)*x(imat,2) + x(j,3)*x(imat,3));
    tmp1=(x(j,1)*dp_tmp - x(imat,1));
    tmp2=(x(j,2)*dp_tmp - x(imat,2));
    tmp3=(x(j,3)*dp_tmp - x(imat,3));
    
    dp = (x(j,1)*x(imat,1) + x(j,2)*x(imat,2) + x(j,3)*x(imat,3));
    
    rd2 = max(0,2*(1-x(imat,1)*x(imat,1).'-x(imat,2)*x(imat,2).'-x(imat,3)*x(imat,3).'));
    rd2v = rd2(:,1);
    
    A(1:fdsize,1:fdsize) = rbf(ep,rd2);
    [LA,UA] = lu(A);

    B(1:fdsize) = (x(j,1)*dp - x(imat,1)).*drbf(ep,rd2v);
    weights = UA\(LA\B);
    
    weights = inv(A)*B;
    weightsDx((j-1)*fdsize+1:j*fdsize) = weights(1:fdsize);

    B(1:fdsize) = (x(j,2)*dp - x(imat,2)).*drbf(ep,rd2v);
    weights = UA\(LA\B);
    weightsDy((j-1)*fdsize+1:j*fdsize) = weights(1:fdsize);

    B(1:fdsize) = (x(j,3)*dp - x(imat,3)).*drbf(ep,rd2v);
    weights = UA\(LA\B);
    weightsDz((j-1)*fdsize+1:j*fdsize) = weights(1:fdsize);

    
    B(1:fdsize) = ep^(2*order)*hyper(ep^2*rd2v,dim,order).*exp(-ep^2*rd2v);
    weights = UA\(LA\B);

    weightsL((j-1)*fdsize+1:j*fdsize) = weights(1:fdsize);
    
end
DPx = sparse(ind_i,ind_j,weightsDx,N,N); 
DPy = sparse(ind_i,ind_j,weightsDy,N,N); 
DPz = sparse(ind_i,ind_j,weightsDz,N,N); 
L = sparse(ind_i,ind_j,weightsL,N,N);

FullDPx=full(DPx)
FullDPy=full(DPy)
FullDPz=full(DPz)
FullL=full(L)


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