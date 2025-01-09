% scatterers phantom with ROIs
% Fig.1 right create an save
close all
clc
clear
addpath(genpath('src'));

crange = [-60,0];
speckleIdx = 5;
parentpath = [pwd '/numerical/'];
saveto = [parentpath 'Test_scatterers/images/'];
scan = linear_scan(linspace(-0.018,0.018,256).', linspace(0.01,0.036+0.01,256).');
path_phantom = 'phantom_9.hdf5';
pht = us_phantom();
pht.read_file(path_phantom);


%p
load([parentpath 'SimulatedData/scatterers/img.mat'], 'img')
figure;imagesc((scan.x_axis)*1e3,(scan.z_axis)*1e3,20*log10(img./max(img(:)))); 
colormap gray; colorbar;
axis equal manual;  
caxis(crange); 
axis off
colorbar off
set(gcf, "Position",[100,100,350,310])
set(gca, 'Position', [0,0,1,1])
pause(0.05)


maskROI = zeros(size(scan.x_matrix));
padROI = 1.8e-3 * 1;
for k=1:size(pht.sca,1)                
    %-- Compute mask inside
    x = pht.sca(k,1);
    z = pht.sca(k,3);                
    %-- Compute mask ROI
    mask = k * ( (scan.x_matrix > (x-padROI)) & ...
                 (scan.x_matrix < (x+padROI)) & ...
                 (scan.z_matrix > (z-padROI)) & ...
                 (scan.z_matrix < (z+padROI)) );
    maskROI = maskROI + mask;                
end
hold on; contour(scan.x_axis*1e3,scan.z_axis*1e3,maskROI,[1 1],'y-');

% save as PDF
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(gcf, 'Renderer', 'opengl');
%print(gcf,[saveto 'pScat'],'-dpdf','-r0')




