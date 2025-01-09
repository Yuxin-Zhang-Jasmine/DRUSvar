% Computing the the observation of Deno 
% by normalizing the meansurements used in DRUS.
% The saved observation By is approximated to be Ix + n.

clear
clc
close all

imageNames = {'simu_reso', 'simu_cont', 'expe_reso', 'expe_cont', 'expe_cross', 'expe_long'};

for i = imageNames
    image = i{1};
    load([pwd '/picmus/Observation/DRUS/' image '.mat']);
    figure; plot(By); title('DRUS')

    scale = max(abs(By));
    By = By ./ scale;
    By = reshape(By, 256, 256);
    By = By';
    By = By(:);
    save([pwd '/picmus/Observation/DENO/' image '.mat'], "By", '-mat')
    figure; plot(By); title('deno') 
end
