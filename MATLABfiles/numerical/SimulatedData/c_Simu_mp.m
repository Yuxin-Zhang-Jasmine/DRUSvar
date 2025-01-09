% load the synthetic echogenicity map (p), 
% add variant multiplicative noise, 
% and save m \odot p (mp)


clear
clc
close all
saveFlag = 1; %!!!!!! [0 or 1]


for type = {'cysts', 'scatterers'}
    type = type{1};
    load([pwd '/numerical/SimulatedData/' type '/img.mat'])
    for i = 1:9
        % gaussian
        x = img.*randn(size(img)) ;
        x = single(x ./ max(abs(x(:))));
        temp = zeros(1,256,256);
        temp(1,:,:) = x ;
    
        x = single(permute(temp, [3,2,1]));
        x = x(:).*4;
    
        if saveFlag
            save([pwd '/numerical/SimulatedData/' type '/' type '_mp_' num2str(i) '.mat'], 'x')
        end
    end

    % view one realization of mp
    img = reshape(single(x ./ max(abs(x))), 256, 256);
    figure; imagesc(20*log10(abs(img))); colormap gray; colorbar;
    axis equal manual; axis([[1,256] [1,256]]); axis off; 
    caxis([-98,0]); 
end