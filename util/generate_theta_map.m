%% Randomly generate surface map of angle theta
clear all;

Ntree = 13; % HTree size in parbrain
num_nvus = sqrt(2^(Ntree-1)); % along one side

N = [num_nvus num_nvus]; % size in pixels of image given by number of NVUs along each side
F = 5;        % frequency-filter width - how big the blobs are
 
% Generate random numbers (normal distribution) then use fft and ifft to
% get blob effect
[X,Y] = ndgrid(1:N(1),1:N(2));
i = min(X-1,N(1)-X+1);
j = min(Y-1,N(2)-Y+1);
H = exp(-.5*(i.^2+j.^2)/F^2);
Z = real(ifft2(H.*fft2(rand(N))));

minZ = min(min(Z)); % For scaling
maxZ = max(max(Z));

theta = pi/(maxZ - minZ) * (Z - minZ);  % Scale so Z is between pi and 0
theta = pi - theta;

% Run once - flips the colorbar
oldcmap = colormap;
colormap( flipud(oldcmap) );

figure(16);
surf(X,Y,theta,'edgecolor','none');
hcb = colorbar;
title(hcb, '\theta')
view([0 90])
xlim([1 num_nvus])
ylim([1 num_nvus])