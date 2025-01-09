% cysts phantom with ROIs
% Fig.1 left create an save
close all
clc
clear
addpath(genpath('src'));

crange = [-60,0];
addpath(genpath('src'));
parentpath = [pwd '/numerical/'];
saveto = [parentpath 'Test_cysts/images/'];
scan = linear_scan(linspace(-0.018,0.018,256).', linspace(0.01,0.036+0.01,256).');
path_phantom = 'phantom_10.hdf5';
pht = us_phantom();
pht.read_file(path_phantom);


% p
load([parentpath 'SimulatedData/cysts/img.mat'], 'img')
fv = 2;  % making the background light gray instead of totally white
img(1,:) = fv;
img(end,:) = fv;
img(:,end) = fv;
img(1,end) = fv;  % should be img(:,1) = fv ?

figure; imagesc((scan.x_axis)*1e3,(scan.z_axis)*1e3,20*log10(img./max(img(:)))); 
colormap gray; colorbar;
axis equal manual; 
caxis(crange); 
axis off
set(gcf, "Position",[100,100,350,310])
set(gca, 'Position', [0,0,1,1])
pause(0.05)


padding = 1.1;
x = scan.x_matrix;
z = scan.z_matrix;  
for k=1:length(pht.occlusionDiameter)
    r = pht.occlusionDiameter(k) / 2;
    rin = r - padding * pht.lateralResolution;
    rout1 = r + padding * pht.lateralResolution;
    rout2 = 1.0*sqrt(rin^2+rout1^2);
    xc = pht.occlusionCenterX(k);
    zc = pht.occlusionCenterZ(k);
    maskOcclusion = ( ((x-xc).^2 + (z-zc).^2) <= r^2);
    maskInside = ( ((x-xc).^2 + (z-zc).^2) <= rin^2);
    maskOutside = ( (((x-xc).^2 + (z-zc).^2) >= rout1^2) & ...
                 (((x-xc).^2 + (z-zc).^2) <= rout2^2) );
%     hold on; contour(scan.x_axis*1e3,scan.z_axis*1e3,maskOcclusion,[1 1],'y-','linewidth',2);
    hold on; contour(scan.x_axis*1e3,scan.z_axis*1e3,maskInside,[1 1],'r-','linewidth',2);
    hold on; contour(scan.x_axis*1e3,scan.z_axis*1e3,maskOutside,[1 1],'g-','linewidth',2);
end


padROIx = pht.RoiPsfTimeX * pht.lateralResolution;
padROIz = pht.RoiPsfTimeZ * pht.axialResolution; 
for k=1:length(pht.RoiCenterX)
    %-- Compute mask inside
    x = pht.RoiCenterX(k);
    z = pht.RoiCenterZ(k);
    %-- Compute mask ROI
    maskROI = k * ( (scan.x_matrix > (x-padROIx(k))) & ...
                 (scan.x_matrix < (x+padROIx(k))) & ...
                 (scan.z_matrix > (z-padROIz(k))) & ...
                 (scan.z_matrix < (z+padROIz(k))) );
    hold on; contour(scan.x_axis*1e3,scan.z_axis*1e3,maskROI,[1 1],'b-','linewidth',2);
end

annotation(gcf,'textbox',...
    [0.295285714285714 0.6 0.0990000000000002 0.0967741935483872],...
    'Color',[0 0 1],...
    'String','ROI',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',14,...
    'FitBoxToText','off');
annotation(gcf,'textbox',...
    [0.181 0.774193548387097 0.0990000000000002 0.0967741935483872],...
    'Color',[1 0 0],...
    'String','in',...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off');
annotation(gcf,'textbox',...
    [0.0695714285714286 0.883870967741936 0.0990000000000002 0.0967741935483871],...
    'Color',[0.466666666666667 0.674509803921569 0.188235294117647],...
    'String',{'out'},...
    'LineStyle','none',...
    'Interpreter','latex',...
    'FontSize',20,...
    'FitBoxToText','off');


% save as PDF
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(gcf, 'Renderer', 'opengl');
%print(gcf,[saveto 'pCysts'],'-dpdf','-r0')



