%% simulate a phantom of scatterers
addpath(genpath('src'));
scan = linear_scan(linspace(-0.018,0.018,256).', linspace(0.01,0.036+0.01,256).');
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

% uncomment to save the echogenicity map 9 times (can be optimized)
for i = 1 : 9
    %save([pwd '/Test_resolution/data/img' num2str(i) '.mat'], 'img')
end
