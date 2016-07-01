% Solves the fifth test case in Cartesian coordinates on the sphere from:
%     Williamson, et. al. "A Standard Test Set for Numerical Approximations to
%     the Shallow Water Equations in Spherical Geometry",  J. Comput. Phys., 
%     102 , 211-224, 1992.

% For details with regard to RBF-FD implemetation of the above test case, see
% Flyer et al., A guide to RBF-generated finite differences for nonlinear transport: '
% Shallow water simulations on a sphere, J. Comput. Phys. 231 (2012) 4078?095

clear all
%fd = 3;    % size of RBF-FD stencil
fd = 7;    % size of RBF-FD stencil

%tend  = 15; % ending time, needs to be in days
tend  = 0.125; % ending time, needs to be in days

order = 4;  % power of Laplacian, L^order
dim = 2 ;   % dimension of stencil, on sphere dim=2

%
% RBF-FD parameters for N = 163842 node set
%
%ep = 14;   % controls width of Gaussian RBF
%dt = 90;   % time step, needs to be in seconds
%gamma = -6.4e-22; % amount of hyperviscosity applied, multiplies Laplacian^order 

%
% RBF-FD parameters for N = 40962 node set
%
% ep = 6.5; % controls width of Gaussian RBF
% dt = 180; % time step, needs to be in seconds
% gamma = -1e-19; % amount of hyperviscosity applied, multiplies Laplacian^order 
%

% RBF-FD parameters for N = 3600 node set
%
% ep = 2; % controls width of Gaussian RBF
% dt = 1000; % time step, needs to be in seconds
% gamma = -2.97e-16; % amount of hyperviscosity applied, multiplies Laplacian^order 

 % RBF-FD parameters for N = 6400 node set
%
 ep = 2.7; % controls width of Gaussian RBF
 dt = 900; % time step, needs to be in seconds
 gamma = -2.98e-17; % amount of hyperviscosity applied, multiplies Laplacian^order 


% Set to dsply=1 if you want to display the time
% Set to plt=1 if you want to plot results at different time-steps.
dsply=1; 
plt=1; 

%fname = load('md/md001.00004');
%fname = load('md/md002.00009');
fname = load('md/md003.00016');
%fname = load('md/md009.00100');

%fname = load('md/md019.00400');
%fname = load('md059.03600');

%fname = load('md079.06400');

%%% Icosehedral
%fname = load('icos40962');
%fname = load('icos163842');
[H,atm] = test_case_5_cart_rk4_fd(fname,ep,fd,order,dim,gamma,dt,tend,dsply,plt);