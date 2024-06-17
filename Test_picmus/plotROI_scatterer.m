function [] = plotROI_scatterer(scan)
%PLOTROI Summary of this function goes here
path_phantom = '~/Documents/MATLAB/01_TMI/src/phantoms/picmus_phantom_3.hdf5';
pht = us_phantom();
pht.read_file(path_phantom);

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

padding = 1;
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

end

