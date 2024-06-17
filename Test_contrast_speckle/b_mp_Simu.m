% cysts phantom
% load the synthetic echogenicity map (p), 
% add variant multiplicative noise, 
% and save m \odot p (mp)

close all
clc
clear

% Signed multiplicative noise
class = 'signed';
for i = 11:19
    load([pwd, '/Test_contrast_speckle/data/img' num2str(i) '.mat'])

    % gaussian
    x = img.*randn(size(img)) ;
    x = single(x ./ max(abs(x(:))));
    temp = zeros(3,256,256);
    for c = 1: 3
    temp(c,:,:) = x ;
    end
    x = single(permute(temp, [3,2,1]));
    x = x(:).*4;

    % uncomment below to save mp
    % save([pwd '/Test_contrast_speckle/data/' class '_m5_' num2str(i) '.mat'], 'x')

end

% view one realization of mp
img = reshape(single(x ./ max(abs(x))), 256, 256, 3);
img = mean(img,3);
figure; imagesc(20*log10(abs(img))); colormap gray; colorbar;
axis equal manual; axis([[1,256] [1,256]]); axis off; 
caxis([-98,0]); 