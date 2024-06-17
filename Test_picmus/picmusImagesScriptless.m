%%% Qualitative results summary (Fig.4)

%% -- phantom-based images loading
close all
clear
clc
addpath(genpath('src'));
DASpath = '../0_otherAlgos/DAS/';
OnePath = '../1_itSensitivity/';
MulPath = '../3_repeatNumSensitivity/';

itLst = [5:5:65, 70:10:100, 120, 140, 160, 200, 250, 333, 500, 1000];
itpicked = 50;
itIdx = find(itLst==itpicked);  % it = '50it';

repeatLst = [3: 13];
repeatIdx = 8;  % corresponds to repeat 10 times (50it x 10)

% parameters
numPhan = 4;
model = {'DAS75', 'DAS1','DENOmean','DRUSmean','DRUSvar'};
numModel = length(model);
scan = linear_scan(linspace(-0.018,0.018,256).', linspace(0.01,0.036+0.01,256).');


% initialization
X = cell(numPhan,numModel);
x_lim = [min(scan.x_matrix(:)) max(scan.x_matrix(:))]*1e3; 
z_lim = [min(scan.z_matrix(:)) max(scan.z_matrix(:))]*1e3;


% load and quantitize each restored image
for i = 1:numModel
    for j = 1:numPhan
        if strcmp(model{i},'DAS1')
            img = us_image;
            img.read_file_hdf5([DASpath, num2str(j), '.hdf5'])
            x = img.data(:,:,1);
        elseif strcmp(model{i},'DAS11')
            img = us_image;
            img.read_file_hdf5([DASpath, num2str(j), '.hdf5'])
            x = img.data(:,:,3);  
        elseif strcmp(model{i},'DAS75')
            img = us_image;
            img.read_file_hdf5([DASpath, num2str(j), '.hdf5'])
            x = img.data(:,:,4);
            
        elseif strcmp(model{i},'DENOone')
            load([OnePath 'Deno/DenoMetricsEnv.mat'])
            x = ImageSet{1, j, itIdx}.data;
        elseif strcmp(model{i},'DRUSone')
            load([OnePath 'DRUS/DRUSMetricsEnv.mat'])
            x = ImageSet{1, j, itIdx}.data;
        elseif strcmp(model{i},'WDRUSone')
            load([OnePath 'WDRUS/WDRUSMetricsEnv.mat'])
            x = ImageSet{1, j, itIdx}.data;
            
        elseif strcmp(model{i},'DENOmean')
            load([MulPath num2str(itpicked) 'it/DenoMetricsmapEnv.mat'])
            x = ImageSetsmap{repeatIdx}{3, j}.data;
        elseif strcmp(model{i},'DRUSmean')
            load([MulPath num2str(itpicked) 'it/DRUSMetricsmapEnv.mat'])
            x = ImageSetsmap{repeatIdx}{3, j}.data;
        elseif strcmp(model{i},'WDRUSMean')
            load([MulPath num2str(itpicked) 'it/WDRUSMetricsmapEnv.mat'])
            x = ImageSetsmap{repeatIdx}{3, j}.data;

        elseif strcmp(model{i},'DenoVar')
            load([MulPath num2str(itpicked) 'it/DenoMetricsmapEnv.mat'])
            x = ImageSetsmap{repeatIdx}{2, j}.data;
        elseif strcmp(model{i},'DRUSvar')
            load([MulPath num2str(itpicked) 'it/DRUSMetricsmapEnv.mat'])
            x = ImageSetsmap{repeatIdx}{2, j}.data;
        elseif strcmp(model{i},'WDRUSVar')
            load([MulPath num2str(itpicked) 'it/WDRUSMetricsmapEnv.mat'])
            x = ImageSetsmap{repeatIdx}{2, j}.data;
        end
        X{j,i} = x;
    end
end

%% -- in-vivo images loading
loadpath = '../5_picmusVivoImg/';

% ----- parameters -----
gammas = [2, 5, 10, 8];   % crossDeno longDeno crossDRUS longDRUS
repeats = [4, 4, 20, 20]; % crossDeno longDeno crossDRUS longDRUS
% diffusion train set
trainset = 'both';  % both (3551 images) is better than cross 
% columns  % it = 50 for all of the samples
models = {'DAS4', 'DAS1', 'DenoMean', 'DRUSMean', 'DRUSVar'};
columnsNum = length(models);
% rows
phantoms = [5,6];  % expe_cross and expe_long
rowsNum = length(phantoms);

totalImagesPerRepeat = 20;  % two phantoms with 10 gammas for each
% ----- Start the calculation -----
% Initialization
Xs = cell(rowsNum, columnsNum);

for rowIdx = 1: rowsNum
    phanIdx = phantoms(rowIdx);
    for columnIdx = 1: columnsNum
        model = models{columnIdx};
        if strcmp(model(1:3), 'DAS')
            img = us_image;
            img.read_file_hdf5([DASpath num2str(phanIdx), '.hdf5'])
            Xs{rowIdx, columnIdx} = img.data(:,:,str2double(model(4))); 
        else
            if strcmp(model(1:4), 'Deno')
                gamma = gammas(1:2);
                repeat = repeats(1:2);
            elseif strcmp(model(1:4), 'DRUS')
                gamma = gammas(3:4);
                repeat = repeats(3:4);
            end

            if strcmp(model(end-2: end), 'One')
                load([loadpath trainset '/' model(1:4) '/' num2str((rowIdx-1)*10+gamma(rowIdx)) '_-1.mat']);
                x1 = x(1,:,:); 
                x2 = x(2,:,:); 
                x3 = x(3,:,:); 
                x = (x1+x2+x3) ./ 3;
                x = squeeze(x);
                Xs{rowIdx, columnIdx} = x;
            else
                tmp = zeros(repeat(rowIdx),256,256);
                for j = 1 : repeat(rowIdx)
                    load([loadpath trainset '/' model(1:4) '/' num2str((j-1)*totalImagesPerRepeat + (rowIdx-1)*10+gamma(rowIdx)) '_-1.mat']);
                    x1 = x(1,:,:); 
                    x2 = x(2,:,:); 
                    x3 = x(3,:,:); 
                    x = (x1+x2+x3) ./ 3;
                    x = squeeze(x);
                    tmp(j,:,:) = x;
                end
                if strcmp(model(end-2: end), 'Var')
                    Xs{rowIdx, columnIdx} = var(tmp,1);
                elseif strcmp(model(end-3: end), 'Mean')
                    Xs{rowIdx, columnIdx} = mean(tmp,1);
                elseif strcmp(model(end-3: end), 'Median')
                    Xs{rowIdx, columnIdx} = median(tmp,1);
                elseif strcmp(model(end-2: end), 'Std')
                    Xs{rowIdx, columnIdx} = std(tmp,1);
                end
            end
        end
    end
end

%% -- plot the reconstructed images with the phantom-ROIs in the 1st col
addpath([pwd '/Test_picmus/']);  % need to use the 'plotROI_xxx.m' functions

figure('Position', [337,372,1200,400+400]); 
h = tiledlayout(4,numModel,'TileSpacing','tight','Padding','none');
for j = [4,3]
    for i = 1:  numModel
        nexttile
        Image_realScale(X{j,i}, '', scan)
        colorbar off; set(gca,'xtick',[]);
        ylabel('Depth [mm]');
        
        % plot the ROIs
        if i == 1
            if j == 4
            plotROI_occlusion(scan)
            else
            plotROI_scatterer(scan)
            end
        end
         
    end
end

for rowIdx = 1:rowsNum
    for columnIdx = 1:columnsNum
        nexttile
        Image_realScale(Xs{rowIdx, columnIdx}, '', scan)
        colorbar off;  
        ylabel('Depth [mm]');
        if rowIdx == 1
            set(gca,'xtick',[]);
        else
            xlabel('Lateral [mm]'); 
        end
    end
end
h.InnerPosition=[0.06 0.055 0.945 0.89]; %[left, bottom, width, height]

%
name = {'75PWs DAS', '1PW DAS', '1PW DENOmean', '1PW DRUSmean', '1PW DRUSvar\bf{(proposed)}'};
for i = 1: numModel
    annstr = name{i}; % annotation text
    annpos = [(i-1)*0.195+0.11 0.95 0.07,0.05]; % annotation position in figure coordinates
    ha = annotation('textbox',annpos,'string',annstr,'Interpreter','latex');
    ha.HorizontalAlignment = 'center';
    ha.BackgroundColor = 'none'; % make the box opaque with some color
    ha.EdgeColor = 'none';
    ha.FontSize = 14;
end

name = {'{\it CL}', '{\it CC}', '{\it ER}', '{\it EC}'};
for i = 1: length(name)
    annstr = (name{i}); % annotation text
    annpos = [ 0.005, (i-1)*0.22+0.18 0.03,0.05]; % annotation position in figure coordinates
    ha = annotation('textbox',annpos,'string',annstr);
    ha.HorizontalAlignment = 'center';
    ha.BackgroundColor = 'none'; % make the box opaque with some color
    ha.EdgeColor = 'none';
    ha.FontSize = 14;
end

%% save as fig / png / PDF
saveto = [pwd '/Test_picmus/images/'];
savefig(gcf, [saveto 'Picmus'])
saveas(gcf, [saveto, 'Picmus'], 'png')

% export as PDF
h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[saveto, 'Picmus'],'-dpdf','-r0')