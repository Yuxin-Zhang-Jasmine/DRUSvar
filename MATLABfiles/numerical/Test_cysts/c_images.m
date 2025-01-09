%%% (cysts) Display the B-mode images in Fig.2 
clc
close all
clear

start=0;
gammaLevels =  [0.02, 0.05, 0.08, 0.11, 0.14, 0.17, 0.2, 0.23, 0.26, 0.29, 0.32, 0.35];
c = length(gammaLevels); % 12
vrange = 60;
repeat = 20;
folder = 'results/sigma1.2';
parentpath = [pwd '/numerical/Test_cysts/'];
filepath = [parentpath folder filesep];
filename = 'm1';
saveto = [parentpath 'images' filesep];

channels = 1;
addpath(genpath('src'));
imgSize = 256;
scan = linear_scan(linspace(-0.018,0.018,imgSize).', linspace(0.01,0.036+0.01,imgSize).');

% load data
xImgs = zeros(repeat, c, imgSize, imgSize);
load([filepath filename filesep '1.mat']); yimgSize = sqrt(length(y_0)/channels);
yImgs = zeros(repeat, c, yimgSize, yimgSize);
for i = 1:c
    for r = 1:repeat
        currentIdx = r  + (start+i-1)*repeat;
        load([filepath filename filesep num2str(currentIdx) '_-1.mat'])
        load([filepath filename filesep num2str(currentIdx) '.mat'])
        yImgs(r, i, :, :) = squeeze(mean(reshape(y_0, yimgSize,yimgSize,channels),3)).';
        xImgs(r, i, :, :) = squeeze(x);
    end
end
yImgs = squeeze(mean(yImgs));
meanImgs = squeeze(mean(xImgs));
varImgs = squeeze(var(xImgs));

% Mean
idx = 1;
to = 3;
for i = [1,to]
    image = us_image(); image.scan = scan; image.number_plane_waves=1;
    tmp = abs(squeeze(meanImgs(i, :, :)));
    image.data = tmp./max(tmp(:));
    image.show(vrange); 
    axis off
    title('')
    colorbar off
    set(gcf, "Position",[100,100,350,350])
    pos = get(gca, 'Position');
    pos(1) = 0; pos(2) = 0;
    pos(3) = 1; pos(4) = 1;
    set(gca, 'Position', [0,0,1,1])
    pause(0.05)
    
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)/3])
    set(gcf, 'Renderer', 'opengl');
    % print(gcf,[saveto 'MeanCysts' num2str(idx)],'-dpdf','-r0')
    idx = idx + 1;
end


% var
idx = 1;
for i = [1,to]
    image = us_image(); image.scan = scan; image.number_plane_waves=1;
    tmp = reshape(abs(squeeze(varImgs(i, :, :))),[imgSize,imgSize]);
    image.data = tmp./max(tmp(:));
    image.show(vrange); 
    axis off
    title('')
    colorbar off
    set(gcf, "Position",[100,100,350,350])
    pos = get(gca, 'Position');
    pos(1) = 0; pos(2) = 0;
    pos(3) = 1; pos(4) = 1;
    set(gca, 'Position', [0,0,1,1])
    pause(0.05)
    
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)/3])
    set(gcf, 'Renderer', 'opengl');
    % print(gcf,[saveto 'VarCysts' num2str(idx)],'-dpdf','-r0')
    idx = idx + 1;
end

% Degraded Images
idx = 1;
for i = [1,to]
    yscan = linear_scan(linspace(-0.018,0.018,yimgSize).', linspace(0.01,0.036+0.01,yimgSize).');
    image = us_image(); image.scan = yscan; image.number_plane_waves=1;
    tmp = abs(squeeze(yImgs(i, :, :)));
    image.data = tmp./max(tmp(:));
    image.show(vrange); 
    axis off
    title('')
    colorbar off
    set(gcf, "Position",[100,100,350,350])
    pos = get(gca, 'Position');
    pos(1) = 0; pos(2) = 0;
    pos(3) = 1; pos(4) = 1;
    set(gca, 'Position', [0,0,1,1])
    pause(0.05)
    
    set(gcf,'Units','Inches');
    pos = get(gcf,'Position');
    set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)/3])
    set(gcf, 'Renderer', 'opengl');
    % print(gcf,[saveto 'yCysts' num2str(idx)],'-dpdf','-r0')
    idx = idx + 1;
end
