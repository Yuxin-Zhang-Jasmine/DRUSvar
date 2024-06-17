%% Compare to the ADMSS in-vivo image (Fig.5)
close all
clear
clc
addpath(genpath('src'));
DASpath = '../0_otherAlgos/DAS/';
OnePath = '../1_itSensitivity/';
MulPath = '../3_repeatNumSensitivity/';

scan = linear_scan(linspace(-0.018,0.018,256).', linspace(0.01,0.036+0.01,256).');
x_lim = [min(scan.x_matrix(:)) max(scan.x_matrix(:))]*1e3; 
z_lim = [min(scan.z_matrix(:)) max(scan.z_matrix(:))]*1e3;
loadpath = '../5_picmusVivoImg/';

gammas = [2, 5, 10, 8];   % crossDeno longDeno crossDRUS longDRUS
repeats = [4, 4, 20, 20]; % crossDeno longDeno crossDRUS longDRUS
% diffusion train set
trainset = 'both';  % both (3551 images) is better than cross 
models = {'DRUSMean', 'DRUSVar'};
columnsNum = length(models);
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
                x = squeeze(mean(x));
                Xs{rowIdx, columnIdx} = x;
            else
                tmp = zeros(repeat(rowIdx),256,256);
                for j = 1 : repeat(rowIdx)
                    load([loadpath trainset '/' model(1:4) '/' num2str((j-1)*totalImagesPerRepeat + (rowIdx-1)*10+gamma(rowIdx)) '_-1.mat']);
                    tmp(j,:,:) = squeeze(mean(x));
                end
                if strcmp(model(end-2: end), 'Var')
                    Xs{rowIdx, columnIdx} = abs(squeeze(var(tmp,1)));
                elseif strcmp(model(end-3: end), 'Mean')
                    Xs{rowIdx, columnIdx} = abs(squeeze(mean(tmp,1)));
                elseif strcmp(model(end-3: end), 'Median')
                    Xs{rowIdx, columnIdx} = abs(squeeze(median(tmp,1)));
                elseif strcmp(model(end-2: end), 'Std')
                    Xs{rowIdx, columnIdx} = abs(squeeze(std(tmp,1)));
                end
            end
        end
    end
end

% DRUSmean + ADMSS 
% AD_param.sigma=0.005;
% AD_param.rho=0.005;
% AD_param.nitmax=19;
% AD_param.n_memory=5;
% AD_param.delta_t=0.3;
AD_param.sigma=4; %4
AD_param.rho=4;
AD_param.nitmax=18; %31
AD_param.n_memory=9;
AD_param.delta_t=0.3;
AD_param.estim='Gamma';
mask = ones(scan.Nz, scan.Nx);

% plot the reconstructed images
addpath([pwd '/Test_picmus/']);

figure('Position', [337,372,440,135]); 
h = tiledlayout(1,2,'TileSpacing','tight','Padding','none');
rowIdx = 1;
for columnIdx = 1:columnsNum
    nexttile
    if columnIdx==1
    Xs{rowIdx, columnIdx} = ADMSS_2D(Xs{rowIdx, columnIdx},mask,AD_param);
    end
    Image_realScale(Xs{rowIdx, columnIdx}, '', scan)
    xlim([-10, 16])
    ylim([10,24])
    colorbar off; axis off
end

h.InnerPosition=[0.0 0.0 1 0.88]; %[left, bottom, width, height]

%% Add titles
name = {'DRUSmean+ADMSS', 'DRUSvar \bf{(proposed)}'};
for i = 1: columnsNum
    annstr = name{i}; % annotation text
    annpos = [(i-1)*0.54 0.97 0.5,0.05]; % annotation position in figure coordinates
    ha = annotation('textbox',annpos,'string',annstr,'Interpreter','latex');
    ha.HorizontalAlignment = 'center';
    ha.BackgroundColor = 'none'; % make the box opaque with some color
    ha.EdgeColor = 'none';
    ha.FontSize = 13;
end

%% Remove all annotations
annotations = findall(gcf,'Type','annotation');
delete(annotations);

%% save as Fig / PNG / PDF
saveto = [pwd filesep 'tryNatural2US/Test_picmus/images/'];
savefig(gcf, [saveto 'edgesCompare'])
saveas(gcf, [saveto, 'edgesCompare'], 'png')

% export as PDF
h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,[saveto, 'edgesCompare'],'-dpdf','-r0')