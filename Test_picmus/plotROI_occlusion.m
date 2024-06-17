function [] = plotROI_occlusion(scan)
%PLOTROI Summary of this function goes here
path_phantom = '~/Documents/MATLAB/01_TMI/src/phantoms/picmus_phantom_4.hdf5';
pht = us_phantom();
pht.read_file(path_phantom);

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

end

