function [DPx,DPy,DPz,L] = rbfmatrix_fd_hyper(x,ep,fdsize,order,dim)
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

for j=1:N
    
    imat = idx(j,:);
    ind_i((j-1)*fdsize+1:j*fdsize) = j;
    ind_j((j-1)*fdsize+1:j*fdsize) = imat;
    
    dp = (x(j,1)*x(imat,1) + x(j,2)*x(imat,2) + x(j,3)*x(imat,3));
    
    rd2 = max(0,2*(1-x(imat,1)*x(imat,1).'-...
        x(imat,2)*x(imat,2).'-x(imat,3)*x(imat,3).'));
    rd2v = rd2(:,1);
    
    A(1:fdsize,1:fdsize) = rbf(ep,rd2);
    [LA,UA] = lu(A);

    B(1:fdsize) = (x(j,1)*dp - x(imat,1)).*drbf(ep,rd2v);
    weights = UA\(LA\B);
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