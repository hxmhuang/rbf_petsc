% RBFInterpolation2D
% Script that performs basic 2D RBF interpolation
% Calls on: CreatePoints, DistanceMatrix, testfunctionsD, PlotSurf,
%           PlotError2D
  % Define the Gaussian RBF and shape parameter
<<<<<<< HEAD:RBFInterpolation2D/orig_matlab_code/RBFInterpolation2D.m
  rbf = @(e,r) exp(-(e*r).^2); ep = 11.1;
=======
  rbf = @(e,r) exp(-(e^2*r)); ep = 11.1;
>>>>>>> 0867a47a9711ce30aff5a9502eb72f985b1d2185:rbf/orig_matlab_code/RBFInterpolation2D.m
  % Create data sites and centers
  m=3;
  n=3;
  dsites = CreatePoints(m*n,2,'u');
  ctrs = dsites;
  % Create evaluation points
  meval=4;
  neval=4;
  epoints = CreatePoints(meval*neval,2,'u');
%  neval = 40; grid = linspace(0,1,neval);
%  [xe,ye] = meshgrid(grid); epoints = [xe(:) ye(:)];

  % Evaluate the test function at the data points
  rhs = testfunctiondD(dsites);
  % Compute distance matrix between the data sites and centers
  DM_data = DistanceMatrix(dsites,ctrs);
  % Compute interpolation matrix
  IM = rbf(ep,DM_data);
  % Compute distance matrix between evaluation points and centers
  DM_eval = DistanceMatrix(epoints,ctrs);
  % Compute evaluation matrix
  EM = rbf(ep,DM_eval);
  % Compute RBF interpolant
  % (evaluation matrix * solution of interpolation system)
  s = EM * (IM\rhs);
  % Compute exact solution, i.e.,
  % evaluate test function on evaluation points
  exact = testfunctiondD(epoints);
  % Compute errors on evaluation grid
  maxerr = norm(s-exact,inf);
  rms_err = norm(s-exact)/neval;
  fprintf('RMS error:     %e\n', rms_err)
  fprintf('Maximum error: %e\n', maxerr)
  % Plot interpolant
  fview = [160,20]; % for Franke's function
  xe = reshape(epoints(:,1),neval,neval);
  ye = reshape(epoints(:,2),neval,neval);
  caption = ['RBF interpolant ',...
      'false colored by absolute error.'];
  PlotSurf(xe,ye,s,neval,exact,maxerr,fview,caption);
  % Plot maximum error
  caption = 'Absolute error for RBF interpolant.';
  PlotError2D(xe,ye,s,exact,maxerr,neval,fview,caption)

