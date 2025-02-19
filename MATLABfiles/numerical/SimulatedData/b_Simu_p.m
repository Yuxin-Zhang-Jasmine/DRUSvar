% create the synthetic echogenicity map (p) and save
% (for cysts and scatterers)

clear
clc
close all
saveFlag = 1; %!!!!!! [0 or 1]


addpath(genpath('src'));
scan = linear_scan(linspace(-0.018,0.018,256).', linspace(0.01,0.036+0.01,256).');

% simulate a phantom of cysts
img = ones(65536, 1);
r  = 3.5/1000;           % Radius of cysts [m] 

xc = [-11, 0, 11, ...
      -11, 0, 11, ...
      -11, 0, 11]./1000; % Places of cysts [m]

zc = [17, 17, 17,  ...
      28, 28, 28,  ...
      39, 39, 39]./1000;


for i = 1 : length(xc)
    inside = ((abs(scan.x-xc(i))).^2 + (abs(scan.z-zc(i))).^2) < r.^2;   
    img = img .* (1-inside);
end

img = reshape(single(img ./ max(abs(img))), 256, 256);
figure; imagesc(20*log10(img)); colormap gray; colorbar;
axis equal manual; axis([[1,256] [1,256]]); axis off; 
caxis([-98,0]); 

if saveFlag
    save([pwd '/numerical/SimulatedData/cysts/img.mat'], 'img')
end



% simulate a phantom of scatterers
img = zeros(65536, 1);
r  = 0.5/1000;                  % Radius of scatters [m] 0.3

xc = [-14, -7, 0, 7, 14, ...
      -14, -7, 0, 7, 14, ...
      -14, -7, 0, 7, 14, ...
      -14, -7, 0, 7, 14, ...
      -14, -7, 0, 7, 14]./1000; % Places of scatters [m]

zc = [14, 14, 14, 14, 14, ...
      21, 21, 21, 21, 21, ...
      28, 28, 28, 28, 28, ...
      35, 35, 35, 35, 35, ...
      42, 42, 42, 42, 42]./1000;


for i = 1 : length(xc)
    inside = ((abs(scan.x-xc(i))).^2 + (abs(scan.z-zc(i))).^2) < r.^2;   
    %inside = (abs(scan.x-xc(i)) < r) .* (abs(scan.z-zc(i)) < r);
    img = img .* (1-inside);
    img = img + inside; 
end

img = reshape(single(img ./ max(abs(img))), 256, 256);
figure; imagesc(20*log10(img)); colormap gray; colorbar;
axis equal manual; axis([[1,256] [1,256]]); axis off; 
caxis([-98,0]); 

if saveFlag
    save([pwd '/numerical/SimulatedData/scatterers/img.mat'], 'img')
end

