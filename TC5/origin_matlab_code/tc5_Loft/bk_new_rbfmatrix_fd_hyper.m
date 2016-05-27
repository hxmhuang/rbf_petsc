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

dp1_tmp = -Px*x'.*MaskMatrix;
dp1=dp1_tmp';
dp2_tmp = -Py*x'.*MaskMatrix;
dp2=dp2_tmp;
dp3_tmp = -Pz*x'.*MaskMatrix;
dp3=dp3_tmp';

%rd2_new= max(0,2*(1-x*x').*MaskMatrix);
rd2_new= max(0,2*(1-x*x.'));
A_new = rbf(ep,rd2_new);

%rd2v_new =drbf(ep,rd2_new);
rd2v_new=-2*ep^2*A_new;

%for j=1:N
for j=1:N
    %rd2v_new = rd2_new(:,j);
    M_tmp1=rd2v_new(:,j);
    %M_tmp2=repmat(M_tmp1,1,N);
    %M_tmp3=M_tmp2';
    B_new1 = dp1(:,j).*M_tmp1;
    B_new2 = dp2(:,j).*M_tmp1;
    B_new3 = dp3(:,j).*M_tmp1;

    K_tmp1= repmat(MaskMatrix(j,:)',1,N);
    K_tmp2= repmat(MaskMatrix(j,:),N,1);
    K_tmp3= K_tmp1.*K_tmp2;
    K_tmp4= K_tmp3|eye(N);
    K_tmp5= A_new.* K_tmp4;
    ExtM=[K_tmp5 MaskMatrix(j,:)'; MaskMatrix(j,:) 0]
    
    RightV=B_new1;
    ExtV=[RightV;0];
    [LA,UA] = lu(ExtM);
    weights1_new = UA\(LA\ExtV);
    DPx_new(j,:)= weights1_new(1:N);
    
    RightV=B_new2;
    ExtV=[RightV;0]
    [LA,UA] = lu(ExtM);
    weights2_new = UA\(LA\ExtV);
    DPy_new(j,:)= weights2_new(1:N);
    
    RightV=B_new3;
    ExtV=[RightV;0]
    [LA,UA] = lu(ExtM);
    weights3_new = UA\(LA\ExtV);
    DPz_new(j,:)= weights3_new(1:N);
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
    P(:,j) = 4*(ep2r2-2*j-d/2+4).*P(:,j-1) - ...
            8*(j-2)*(2*j+d-6)*P(:,j-2);
end
p = P(:,k+1);