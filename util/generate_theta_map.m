%% Generate surface maps of angle theta for use in parbrain
% To use: 
% 1. Generate desired theta matrix
% 2. Save theta variable as csv: csvwrite('filename', theta)
% 3. Use in the C code as a command line argument - see C code for details

clear all;

%% For use in all maps
Ntree = 13;                     % Number of levels in the H-tree size in parbrain
num_nvus = sqrt(2^(Ntree-1));   % number of NVUs along one side
lim = num_nvus/2 - 0.5;         % For plotting

% Generate X,Y coordinate grid
[X,Y] = meshgrid(-lim:lim,-lim:lim);

% Run once - flips the colorbar around
% oldcmap = colormap;
% colormap( flipud(oldcmap) );

%% Random surface generator
% For modelling and simulative purposes random rough surfaces with Gaussian statistics 
% can be generated using a method outlined by Garcia and Stoll (1984), where an 
% uncorrelated distribution of surface points using a random number generator 
% (i.e. white noise) is convolved with a Gaussian filter to achieve correlation. 
% This convolution is most efficiently performed using the discrete 
% Fast Fourier Transform (FFT) algorithm.

% Generates a square 2-dimensional random rough surface Z(X,Y) with num_nvus x num_nvus 
% surface points. The surface has a Gaussian height distribution and 
% exponential autocovariance functions which describes the covariance (correlation) 
% of the surface with translationally shifted versions, where F is the 
% correlation length, i.e. the typical distance between two similar features 
% (hills or valleys). 

F = 5;        % correlation length - how big the blobs are
 
% The exponential covariance function (Gaussian filter)
H = exp(-.5*(X.^2+Y.^2)/F^2);

% Random grid of numbers
z = rand(num_nvus);

% Convolve the random numbers with the covariance function
Z = ifft2(fft2(H).*fft2(z));

% Scale so theta is between pi and 0
minZ = min(min(Z)); 
maxZ = max(max(Z));
theta = pi*((maxZ - Z)./(maxZ-minZ));  

% Plot curvature 
r = 20/(2*pi);  % Minor radius
n = 4;        % Ratio n=R/r
R = n*r;

a = r*sqrt(n^2 - 1);
eta = atanh(a/(n*r));

curve = 1/(r^2) - (n*(n-cos(theta)) ) ./ (a^2);

figure(17);
surf(X,Y,curve,'edgecolor','none');
hcb = colorbar;
title(hcb, {'$\Gamma(\tilde{\theta})$'},'Interpreter','latex')
view([0 90])
xlim([-lim lim])
ylim([-lim lim])
caxis([-0.033 0.02])
xlabel('x');
ylabel('y');

figure(170);
surf(X,Y,theta,'edgecolor','none');
hcb = colorbar;
title(hcb, {'$\tilde{\theta}$'},'Interpreter','latex')
view([0 90])
xlim([-lim lim])
ylim([-lim lim])
caxis([0 pi])
xlabel('x');
ylabel('y');

%% Theta mapping with cos(x) - tube in the middle of the tissue slice

% function that varies sinusoidally with the x coordinate and is scaled for
% theta between 0 and pi
theta2 = (pi/2)*(1+cos(X/6));

% Plot curvature 
r = 20/(2*pi);  % Minor radius
n = 4;        % Ratio n=R/r
R = n*r;

a = r*sqrt(n^2 - 1);
eta = atanh(a/(n*r));

curve2 = 1/(r^2) - (n*(n-cos(theta2)) ) ./ (a^2);

figure(18);
surf(X,Y,curve2,'edgecolor','none');
hcb = colorbar;
title(hcb, {'$\Gamma(\tilde{\theta})$'},'Interpreter','latex')
view([0 90])
xlim([-lim lim])
ylim([-lim lim])
caxis([-0.033 0.02])
xlabel('x');
ylabel('y');

figure(180);
surf(X,Y,theta2,'edgecolor','none');
hcb = colorbar;
title(hcb, {'$\tilde{\theta}$'},'Interpreter','latex')
view([0 90])
xlim([-lim lim])
ylim([-lim lim])
caxis([0 pi])
xlabel('x');
ylabel('y');


%% Mapping with positive curvature sections - 2 ellipsoids

% Ellipsoid parameters
xr = 30; yr = 10; zr = 1;   % Scale 
x0 = -lim; x1=lim;          % Centers of ellipsoids
y0 = -15;                   % Move up or down

% Two ellipsoids centered at x0,y0 and x1,y0 
Z3 = real(sqrt(  xr^2 * (yr^2 - (y0 - Y).^2) - yr^2 * (x0 - X).^2  )./(xr*yr)) + real(sqrt(  xr^2 * (yr^2 - (y0 - Y).^2) - yr^2 * (x1 - X).^2  )./(xr*yr) );

% Scale so theta is between pi and 0
theta3 = pi*(1 - Z3);

% Set the remainder of the tissue slice to almost flat (slightly negative curvature)
theta3 = min(pi/2, theta3);

% Plot curvature 
r = 20/(2*pi);  % Minor radius
n = 4;        % Ratio n=R/r
R = n*r;

a = r*sqrt(n^2 - 1);
eta = atanh(a/(n*r));

curve3 = 1/(r^2) - (n*(n-cos(theta3)) ) ./ (a^2);

figure(19);
surf(X,Y,curve3,'edgecolor','none');
hcb = colorbar;
title(hcb, {'$\Gamma(\tilde{\theta})$'},'Interpreter','latex')
view([0 90])
xlim([-lim lim])
ylim([-lim lim])
caxis([-0.033 0.02])
xlabel('x');
ylabel('y');

figure(190);
surf(X,Y,theta3,'edgecolor','none');
hcb = colorbar;
title(hcb, {'$\tilde{\theta}$'},'Interpreter','latex')
view([0 90])
xlim([-lim lim])
ylim([-lim lim])
caxis([0 pi])
xlabel('x');
ylabel('y');
