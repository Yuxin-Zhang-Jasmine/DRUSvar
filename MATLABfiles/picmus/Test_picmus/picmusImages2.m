%%% Qualitative results summary (Fig.4)

close all
clear
clc
addpath(genpath('src'));
DASpath = [pwd '/picmus/DAS/'];
DENOvitropath = [pwd, '/picmus/Test_picmus/results/DENOvitro/'];
DRUSvitropath = [pwd, '/picmus/Test_picmus/results/DRUSvitro/'];
DENOvivopath = [pwd, '/picmus/Test_picmus/results/DENOvivo/'];
DRUSvivopath = [pwd, '/picmus/Test_picmus/results/DRUSvivo/'];

repeat = 10;
phantoms = {'simu_reso', 'simu_cont', 'expe_reso', 'expe_cont', 'expe_cross', 'expe_long'};
phanSelected = [4,3,5,6];  % show 'expe_cont','expe_reso','expe_cross','expe_long'
numRows = length(phanSelected);  
methods = {'DAS75', 'DAS1','DENOmean','DRUSmean','DRUSvar'}; % different approaches
numMethods = length(methods);
scan = linear_scan(linspace(-0.018,0.018,256).', linspace(0.01,0.036+0.01,256).');

% initialization
Xs = cell(numRows,numMethods);
rowId = 0;

% load each restored image
for phanId = phanSelected 
    rowId = rowId + 1;
        % DAS75 & DAS1
        img = us_image;
        img.read_file_hdf5([DASpath, num2str(phanId), '.hdf5'])
        Xs{rowId,1} = img.data(:,:,4);     
        Xs{rowId,2} = img.data(:,:,1);
  
        % DENOmean
        temp = zeros(repeat, 256, 256);
        if phanId < 5 % vitro
            for c = 1:repeat
                load([DENOvitropath, num2str(phanId + (c-1)*4), '_-1.mat']);
                temp(c,:,:) = squeeze(x);
            end
        else          % vivo
            for c = 1:repeat
                load([DENOvivopath, num2str(phanId-4 + (c-1)*2), '_-1.mat']);
                temp(c,:,:) = squeeze(x);
            end
        end
        Xs{rowId,3} = mean(temp); 
        
        % DRUSmean & DRUSvar
        temp = zeros(repeat, 256, 256);
        if phanId < 5 % vitro
            for c = 1:repeat
                load([DRUSvitropath, num2str(phanId + (c-1)*4), '_-1.mat']);
                temp(c,:,:) = squeeze(x);
            end
        else          % vivo
            for c = 1:repeat
                load([DRUSvivopath, num2str(phanId-4 + (c-1)*2), '_-1.mat']);
                temp(c,:,:) = squeeze(x);
            end
        end            
        Xs{rowId,4} = squeeze(mean(temp)); 
        Xs{rowId,5} = squeeze(var(temp));
end


% -- plot the reconstructed images with the phantom-ROIs in the 1st col
addpath([pwd '/picmus/Test_picmus/']);  % need to use 'plotROI_xxx.m' functions

figure('Position', [337,372,1200,400+400]); 
h = tiledlayout(numRows,numMethods,'TileSpacing','tight','Padding','none');
for rowId = 1: numRows
    for i = 1: numMethods
        nexttile
        Image_realScale(Xs{rowId,i}, '', scan)
        colorbar off;
        ylabel('Depth [mm]');
        if rowId == numRows
            xlabel('Lateral [mm]'); 
        else
            set(gca,'xtick',[]);
        end
        % plot the ROIs
        if i == 1
            if rowId == 1
            plotROI_occlusion(scan)
            elseif rowId == 2
            plotROI_scatterer(scan)
            end
        end

    end
end

h.InnerPosition=[0.06 0.055 0.945 0.89]; %[left, bottom, width, height]
name = {'75PWs DAS', '1PW DAS', '1PW DENOmean', '1PW DRUSmean', '1PW DRUSvar\bf{(proposed)}'};
for i = 1: numMethods
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
saveto = [pwd 'numerical/picmus/Test_picmus/images/'];
%savefig(gcf, [saveto 'Picmus'])
%saveas(gcf, [saveto, 'Picmus'], 'png')

% export as PDF
h = gcf;
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(gcf, 'Renderer', 'opengl');
%print(h,[saveto, 'Picmus'],'-dpdf','-r0')