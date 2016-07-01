%
% Set up for the Williamson test case 5.
%
function atm = setupT5(nfile)

nodes = nfile;
atm.pts.x = nodes(:,1); 
atm.pts.y = nodes(:,2); 
atm.pts.z = nodes(:,3);
atm.pts.nd = size(atm.pts.x,1);
[atm.pts.la,atm.pts.th,r] = cart2sph(nodes(:,1),nodes(:,2),nodes(:,3));

% Variables for projecting an arbitrary Cartesian vector onto the surface
% of the sphere.
x2 = nodes(:,1).^2; xy = nodes(:,1).*nodes(:,2);
y2 = nodes(:,2).^2; xz = nodes(:,1).*nodes(:,3);
z2 = nodes(:,3).^2; yz = nodes(:,2).*nodes(:,3);
atm.pts.p_u = [1-x2  -xy   -xz];
atm.pts.p_v = [-xy  1-y2   -yz];
atm.pts.p_w = [-xz   -yz  1-z2];

atm.alpha = 0;         % Angle of rotation measured from the equator.
atm.u0 = 20;           % Speed of rotation in meters/second
atm.a = 6.37122e6;     % Mean radius of the earth (meters).
atm.omega = 7.292e-5;  % Rotation rate of the earth (1/seconds).
atm.g = 9.80616;       % Gravitational constant (m/s^2).
atm.gh0 = atm.g*5960;  % Initial condition for the geopotential field (m^2/s^2).
atm.f = 2*atm.omega*(-atm.pts.x*sin(atm.alpha) + atm.pts.z*cos(atm.alpha)); % Coriolis force.

% Parameters for the mountain:
atm.lam_c = -pi/2;
atm.thm_c = pi/6;
atm.mR = pi/9;
atm.hm0 = 2000;

% Compute the profile of the mountain (multiplied by gravity)
atm.ghm = zeros(atm.pts.nd,1);
r2 = (atm.pts.la-atm.lam_c).^2 + (atm.pts.th-atm.thm_c).^2;
id = r2 < atm.mR^2;
atm.ghm(id) = atm.g*atm.hm0*(1-sqrt(r2(id))/atm.mR);
