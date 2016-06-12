function [H,atm] = test_case_5_cart_rk4_fd(nfile,ep,fdsize,order,dim,gamma,dt,tend,dsply,plt)

% Initialize the constants for the Williamson test case 5.
atm = setupT5(nfile);
x = atm.pts.x; 
y = atm.pts.y; 
z = atm.pts.z;


% Calculate the RBF-FD differentiation matrices
[DPx,DPy,DPz,L] = rbfmatrix_fd_hyper([x y z],ep,fdsize,order,dim);
L = gamma*L;


% Initial condition.
[uc,gh] = computeInitialCondition(atm);
% H contains the velocity components in x,y,z and the geopotential without
% the mean offset gh0.
H = [uc gh];


% Compute the projected gradient of the mountain for test case 5
gradghm = zeros(atm.pts.nd,3);
gradghm(:,1) = DPx*atm.ghm/atm.a;
gradghm(:,2) = DPy*atm.ghm/atm.a;
gradghm(:,3) = DPz*atm.ghm/atm.a;


% For plotting results during time-stepping
[lam,theta] = cart2sph(x,y,z);  
tri_lt = delaunay(lam,theta);  
dsply_rate = 10;       % Number of time-steps to wait before plotting a result.

timeday = zeros((3600/dt)*24*tend,1); % create space in memory for the time in days

% Time-step using RK4
for nt=1:tend*24*3600
   K = H;
   d1 = dt*evalCartRhs_fd(K,DPx,DPy,DPz,L,atm,gradghm);
   K = H + 0.5*d1;
   d2 = dt*evalCartRhs_fd(K,DPx,DPy,DPz,L,atm,gradghm);
   K = H + 0.5*d2;
   d3 = dt*evalCartRhs_fd(K,DPx,DPy,DPz,L,atm,gradghm);
   K = H + d3;
   d4 = dt*evalCartRhs_fd(K,DPx,DPy,DPz,L,atm,gradghm);
   H = H + 1/6*(d1 + 2*d2 + 2*d3 + d4);

   % Record the current time in days.
   timeday(nt,1) = (nt*dt)/(3600*24);

   % Print out time in days and plot after every dsply_rate time steps.
   if mod(nt,dsply_rate) == 0 || nt == 1
      if dsply == 1
         fprintf('time\n   ');
         fprintf('%1.3e\n\n  ',timeday(nt));
         if plt == 1
            tricontour(tri_lt,lam,theta,(H(:,4)+atm.gh0)/atm.g,linspace(5e3,6e3,11));
            axis([-pi-0.01 pi+0.01 -pi/2-0.01 pi/2+0.01]);
            xlabel('\lambda');
            ylabel('\theta');
            pause(0.01);
         end
      end
   end

end

