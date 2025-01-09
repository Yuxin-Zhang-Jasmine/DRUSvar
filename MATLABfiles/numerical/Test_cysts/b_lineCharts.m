%%% (cysts) Display the line charts in Fig.2 
%% evaluate and save the metrics
clc
close all
clear

flag_display = 0;
flag_avg = 0;
channels = 1;
addpath(genpath('src'));
path_phantom = 'phantom_10.hdf5';
imgSize = 256;
scan = linear_scan(linspace(-0.018,0.018,imgSize).', linspace(0.01,0.036+0.01,imgSize).');
c = 12;  % 12 levels of additive noise
repeat = 20;
REPEAT = 20;
m = [1, 2, 3, 4, 5, 6, 7, 8];
filepath = [pwd '/numerical/Test_cysts/results/sigma1.2/'];
metrics = cell(length(m),4,7);

for p = 1:length(m)
    filename = ['m' num2str(m(p))];
    % load data
    xImgs = zeros(repeat, c, imgSize, imgSize);
    load([filepath filename filesep '1.mat']); yimgSize = sqrt(length(y_0)/channels);
    yImgs = zeros(repeat, c, yimgSize, yimgSize);
    for i = 1:c
        for r = 1:repeat
            currentIdx = r  + (i-1)*REPEAT;
            load([filepath filename filesep num2str(currentIdx) '_-1.mat'])
            load([filepath filename filesep num2str(currentIdx) '.mat'])
            yImgs(r, i, :, :) = squeeze(mean(reshape(y_0, yimgSize,yimgSize,channels),3)).';
            xImgs(r, i, :, :) = squeeze(x);
        end
    end
    yImgs = squeeze(mean(yImgs));
    meanImgs = squeeze(mean(xImgs));
    varImgs = squeeze(var(xImgs)); 
    
    
    
    % Degraded Images
    for i = 1: c
    yscan = linear_scan(linspace(-0.018,0.018,yimgSize).', linspace(0.01,0.036+0.01,yimgSize).');
    image = us_image(); image.scan = yscan; image.number_plane_waves=1;
    tmp = abs(squeeze(yImgs(i, :, :)));
    image.data = tmp./max(tmp(:)); 
    [yFWHMA(i,:), yFWHML(i,:), yCNR(i,:), ygCNR(i,:), ySNR(i,:), yKS(i,:)] = evaluation(path_phantom, image, flag_display, flag_avg); 
    % image.postenv = 1;
    % [yFWHMA(i,:), yFWHML(i,:), yCNR(i,:), ygCNR(i,:), ySNR(i,:), yKS(i,:)] = evaluation(path_phantom, image, flag_display, flag_avg); 
    % image.postenv = 2;
    % [~, yFWHML(i,:), yCNR(i,:), ygCNR(i,:), ySNR(i,:), yKS(i,:)] = evaluation(path_phantom, image, flag_display, flag_avg); 
    end
    
    % Mean
    for i = 1: c
    image = us_image(); image.scan = scan; image.number_plane_waves=1;
    tmp = abs(squeeze(meanImgs(i, :, :)));
    image.data = tmp./max(tmp(:)); 
    [MeanFWHMA(i,:), MeanFWHML(i,:), MeanCNR(i,:), MeangCNR(i,:), MeanSNR(i,:), MeanKS(i,:)] = evaluation(path_phantom, image, flag_display, flag_avg); 
    % image.postenv = 1;
    % [MeanFWHMA(i,:), MeanFWHML(i,:), MeanCNR(i,:), MeangCNR(i,:), MeanSNR(i,:), MeanKS(i,:)] = evaluation(path_phantom, image, flag_display, flag_avg); 
    % image.postenv = 2;
    % [~, MeanFWHML(i,:), MeanCNR(i,:), MeangCNR(i,:), MeanSNR(i,:), MeanKS(i,:)] = evaluation(path_phantom, image, flag_display, flag_avg); 
    end
    
    % var
    for i = 1: c
    image = us_image(); image.scan = scan; image.number_plane_waves=1;
    tmp = reshape(abs(squeeze(varImgs(i, :, :))),[imgSize,imgSize]);
    image.data = tmp./max(tmp(:)); 
    [VarFWHMA(i,:), VarFWHML(i,:), VarCNR(i,:), VargCNR(i,:), VarSNR(i,:), VarKS(i,:)] = evaluation(path_phantom, image, flag_display, flag_avg);
    % image.postenv = 1;
    % [VarFWHMA(i,:), VarFWHML(i,:), VarCNR(i,:), VargCNR(i,:), VarSNR(i,:), VarKS(i,:)] = evaluation(path_phantom, image, flag_display, flag_avg);
    % image.postenv = 2;
    % [~, VarFWHML(i,:), VarCNR(i,:), VargCNR(i,:), VarSNR(i,:), VarKS(i,:)] = evaluation(path_phantom, image, flag_display, flag_avg);
    end
    
    metrics{p,1,1} = MeanFWHMA;
    metrics{p,1,2} = MeanFWHML;
    metrics{p,1,3} = MeanCNR;
    metrics{p,1,4} = MeangCNR;
    metrics{p,1,5} = MeanSNR;
    metrics{p,1,6} = MeanKS;
    metrics{p,1,7} = "Mean";
    
    metrics{p,2,1} = VarFWHMA;
    metrics{p,2,2} = VarFWHML;
    metrics{p,2,3} = VarCNR;
    metrics{p,2,4} = VargCNR;
    metrics{p,2,5} = VarSNR;
    metrics{p,2,6} = VarKS;
    metrics{p,2,7} = "Var";
    
    metrics{p,3,1} = "FWHM_{axial} [mm]";
    metrics{p,3,2} = "FWHM_{lateral} [mm]";
    metrics{p,3,3} = "CNR [dB]";
    metrics{p,3,4} = "gCNR";
    metrics{p,3,5} = "SNR";
    metrics{p,3,6} = "KS";
    
    
    metrics{p,4,1} = yFWHMA;
    metrics{p,4,2} = yFWHML;
    metrics{p,4,3} = yCNR;
    metrics{p,4,4} = ygCNR;
    metrics{p,4,5} = ySNR;
    metrics{p,4,6} = yKS;
    metrics{p,4,7} = "y";
end
% save([pwd '/numerical/Test_cysts/results/metrics.mat'], 'metrics', 'c'); 

%% load the metrics and plot
clc
load([pwd '/numerical/Test_cysts/results/metrics.mat'])
flag_y=1;
type = 1;  fontsize=19;

[ml, ~, ~] = size(metrics);
% 9 cysts in one multiplicative noise realization; 
% 4 ROIs for evaluating SNR
% ml=7 realizations of multiplicative noise (2 are discarded)
% 3 types of images (the degraded one , DRUSmean, DRUSvar)
% c = 12 levels of additive noise
CNR = zeros(9*ml,3,c);
gCNR = zeros(9*ml,3,c);
SNR = zeros(4*ml,3,c);
for p = 1:ml
    for q = 1: 9
    CNR(9*(p-1)+q, 1,:) = metrics{p,1,3}(:,q);
    CNR(9*(p-1)+q, 2,:) = metrics{p,2,3}(:,q);
    CNR(9*(p-1)+q, 3,:) = metrics{p,4,3}(:,q);
    gCNR(9*(p-1)+q, 1,:) = metrics{p,1,4}(:,q);
    gCNR(9*(p-1)+q, 2,:) = metrics{p,2,4}(:,q);
    gCNR(9*(p-1)+q, 3,:) = metrics{p,4,4}(:,q);
    end
    for q = 1: 4
    SNR(4*(p-1)+q, 1,:) = metrics{p,1,5}(:,q);
    SNR(4*(p-1)+q, 2,:) = metrics{p,2,5}(:,q);
    SNR(4*(p-1)+q, 3,:) = metrics{p,4,5}(:,q);
    end    
end
CNR(isnan(CNR)) = 0;
gCNR(isnan(gCNR)) = 0;
SNR(isnan(SNR)) = 0;

meanCNR = squeeze(mean(CNR)); meangCNR = squeeze(mean(gCNR)); meanSNR = squeeze(mean(SNR));
stdCNR = squeeze(std(CNR)); stdgCNR = squeeze(std(gCNR));   stdSNR = squeeze(std(SNR));

% compute the bounds of the std
CNRlowerBound1 = flip((meanCNR - stdCNR),2);
CNRupperBound1 = meanCNR + stdCNR;
CNRlowerBound2 = flip((meanCNR - 2 * stdCNR),2);
CNRupperBound2 = meanCNR + 2 * stdCNR;
CNRlowerBound3 = flip((meanCNR - 3 * stdCNR),2);
CNRupperBound3 = meanCNR + 3 * stdCNR;

gCNRlowerBound1 = flip((meangCNR - stdgCNR),2);
gCNRupperBound1 = meangCNR + stdgCNR;
gCNRlowerBound2 = flip((meangCNR - 2 * stdgCNR),2);
gCNRupperBound2 = meangCNR + 2 * stdgCNR;
gCNRlowerBound3 = flip((meangCNR - 3 * stdgCNR),2);
gCNRupperBound3 = meangCNR + 3 * stdgCNR;

SNRlowerBound1 = flip((meanSNR - stdSNR),2);
SNRupperBound1 = meanSNR + stdSNR;
SNRlowerBound2 = flip((meanSNR - 2 * stdSNR),2);
SNRupperBound2 = meanSNR + 2 * stdSNR;
SNRlowerBound3 = flip((meanSNR - 3 * stdSNR),2);
SNRupperBound3 = meanSNR + 3 * stdSNR;

% prepare the x-axis coordinates for drawing
x = [0.02, 0.05, 0.08, 0.11, 0.14, 0.17, 0.2, 0.23, 0.26, 0.29, 0.32, 0.35];
xconf = [x flipud(x')'];

CNRyconf1 = [CNRupperBound1 CNRlowerBound1];
CNRyconf2 = [CNRupperBound2 CNRlowerBound2];
CNRyconf3 = [CNRupperBound3 CNRlowerBound3];

gCNRyconf1 = [gCNRupperBound1 gCNRlowerBound1];
gCNRyconf2 = [gCNRupperBound2 gCNRlowerBound2];
gCNRyconf3 = [gCNRupperBound3 gCNRlowerBound3];

SNRyconf1 = [SNRupperBound1 SNRlowerBound1];
SNRyconf2 = [SNRupperBound2 SNRlowerBound2];
SNRyconf3 = [SNRupperBound3 SNRlowerBound3];

% chart of CNR
figure('Position', [476 356 700 370] );hold on;
switch type    
    case 1        
        fill(xconf, CNRyconf1(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h1=plot(x, meanCNR(1,:), 'Color', [0 72 186]/255, 'LineStyle', '-', 'LineWidth', 2);  
        plot(x,meanCNR(1,:),'o','MarkerSize',8,'MarkerFaceColor',[0 72 186]/255,'MarkerEdgeColor',[0 72 186]/255,'HandleVisibility', 'off')        

        fill(xconf, CNRyconf1(2,:), [219 112 147]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h2=plot(x, meanCNR(2,:), 'Color', [219 112 147]/255, 'LineStyle', '-', 'LineWidth', 2); 
        plot(x,meanCNR(2,:),'o','MarkerSize',8,'MarkerFaceColor',[219 112 147]/255,'MarkerEdgeColor',[219 112 147]/255,'HandleVisibility', 'off')        
        if flag_y == 1
        fill(xconf, CNRyconf1(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h3=plot(x, meanCNR(3,:), 'Color', [119 172 48]/255, 'LineStyle', '-', 'LineWidth', 2);  
        plot(x,meanCNR(3,:),'o','MarkerSize',8,'MarkerFaceColor',[119 172 48]/255,'MarkerEdgeColor',[119 172 48]/255,'HandleVisibility', 'off')        
        end
    case 2        
        fHdl(1) = fill(xconf, CNRyconf1(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, CNRyconf2(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, CNRyconf3(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h1=plot(x, meanCNR(1,:), 'Color', [0 72 186]/255, 'LineStyle', '-', 'LineWidth', 2);    
        fHdl(1) = fill(xconf, CNRyconf1(2,:), [219 112 147]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, CNRyconf2(2,:), [219 112 147]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, CNRyconf3(2,:), [219 112 147]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h2=plot(x, meanCNR(2,:), 'Color', [219 112 147]/255, 'LineStyle', '-', 'LineWidth', 2);  
        fHdl(1) = fill(xconf, CNRyconf1(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, CNRyconf2(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, CNRyconf3(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h3=plot(x, meanCNR(3,:), 'Color', [119 172 48]/255, 'LineStyle', '-', 'LineWidth', 2);            
end
ax = gca;
ax.Box = 'off';
ax.FontName = 'Times New Roman';
ax.GridLineStyle = '-.';
ax.GridColor = 'k';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.LineWidth = 1;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.TickDir = 'in';
ax.FontSize = fontsize-3;
ylabel('CNR [dB]', 'Interpreter','latex',FontSize=fontsize+2)
xlabel('$std(\mathbf{n})$', 'Interpreter','latex',FontSize=fontsize+2)
legend([h2, h1, h3],["DRUSvar", "DRUSmean", "$\mathbf{By}$"], ...
    "Interpreter","latex", "Orientation","vertical", "Location","best", ...
    fontsize=fontsize);
legend boxoff
xlim([0.018,max(x)]);


% chart of gCNR
figure('Position', [476 356 700 370] );hold on;
switch type    
    case 1        
        fill(xconf, gCNRyconf1(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h1=plot(x, meangCNR(1,:), 'Color', [0 72 186]/255, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName','Mean');  
        plot(x,meangCNR(1,:),'o','MarkerSize',8,'MarkerFaceColor',[0 72 186]/255,'MarkerEdgeColor',[0 72 186]/255,'HandleVisibility', 'off')        

        fill(xconf, gCNRyconf1(2,:), [219 112 147]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h2=plot(x, meangCNR(2,:), 'Color', [219 112 147]/255, 'LineStyle', '-', 'LineWidth', 2,'DisplayName','Var'); 
        plot(x,meangCNR(2,:),'o','MarkerSize',8,'MarkerFaceColor',[219 112 147]/255,'MarkerEdgeColor',[219 112 147]/255,'HandleVisibility', 'off')        
        if flag_y == 1
        fill(xconf, gCNRyconf1(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h3=plot(x, meangCNR(3,:), 'Color', [119 172 48]/255, 'LineStyle', '-', 'LineWidth', 2,'DisplayName','y');  
        plot(x,meangCNR(3,:),'o','MarkerSize',8,'MarkerFaceColor',[119 172 48]/255,'MarkerEdgeColor',[119 172 48]/255,'HandleVisibility', 'off')        
        end
    case 2        
        fHdl(1) = fill(xconf, gCNRyconf1(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, gCNRyconf2(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, gCNRyconf3(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h1=plot(x, meangCNR(1,:), 'Color', [0 72 186]/255, 'LineStyle', '-', 'LineWidth', 2);    
        fHdl(1) = fill(xconf, gCNRyconf1(2,:), [219 112 147]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, gCNRyconf2(2,:), [219 112 147]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, gCNRyconf3(2,:), [219 112 147]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h2=plot(x, meangCNR(2,:), 'Color', [219 112 147]/255, 'LineStyle', '-', 'LineWidth', 2);  
        fHdl(1) = fill(xconf, gCNRyconf1(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, gCNRyconf2(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, gCNRyconf3(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h3=plot(x, meangCNR(3,:), 'Color', [119 172 48]/255, 'LineStyle', '-', 'LineWidth', 2);            
end
ax = gca;
ax.Box = 'off';
ax.FontName = 'Times New Roman';
ax.GridLineStyle = '-.';
ax.GridColor = 'k';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.LineWidth = 1;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.TickDir = 'in';
ax.FontSize = fontsize-3;
ylabel('gCNR', 'Interpreter','latex',FontSize=fontsize+2)
xlabel('$std(\mathbf{n})$', 'Interpreter','latex',FontSize=fontsize+2)
legend([h2, h1, h3],["DRUSvar", "DRUSmean", "$\mathbf{By}$"], ...
    "Interpreter","latex", "Orientation","vertical", "Location","best", ...
    fontsize=fontsize);
legend boxoff
xlim([0.018,max(x)]);
ylim([0,1.05])


% chart of SNR
figure('Position', [476 356 700 370]);hold on;
switch type    
    case 1        
        fill(xconf, SNRyconf1(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h1=plot(x, meanSNR(1,:), 'Color', [0 72 186]/255, 'LineStyle', '-', 'LineWidth', 2);  
        plot(x,meanSNR(1,:),'o','MarkerSize',8,'MarkerFaceColor',[0 72 186]/255,'MarkerEdgeColor',[0 72 186]/255,'HandleVisibility', 'off')        

        fill(xconf, SNRyconf1(2,:), [219 112 147]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h2=plot(x, meanSNR(2,:), 'Color', [219 112 147]/255, 'LineStyle', '-', 'LineWidth', 2); 
        plot(x,meanSNR(2,:),'o','MarkerSize',8,'MarkerFaceColor',[219 112 147]/255,'MarkerEdgeColor',[219 112 147]/255,'HandleVisibility', 'off')        
        if flag_y == 1
        fill(xconf, SNRyconf1(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h3=plot(x, meanSNR(3,:), 'Color', [119 172 48]/255, 'LineStyle', '-', 'LineWidth', 2);  
        plot(x,meanSNR(3,:),'o','MarkerSize',8,'MarkerFaceColor',[119 172 48]/255,'MarkerEdgeColor',[119 172 48]/255,'HandleVisibility', 'off')        
        end
    case 2        
        fHdl(1) = fill(xconf, SNRyconf1(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, SNRyconf2(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, SNRyconf3(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h1=plot(x, meanSNR(1,:), 'Color', [0 72 186]/255, 'LineStyle', '-', 'LineWidth', 2);    
        fHdl(1) = fill(xconf, SNRyconf1(2,:), [219 112 147]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, SNRyconf2(2,:), [219 112 147]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, SNRyconf3(2,:), [219 112 147]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h2=plot(x, meanSNR(2,:), 'Color', [219 112 147]/255, 'LineStyle', '-', 'LineWidth', 2);  
        fHdl(1) = fill(xconf, SNRyconf1(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, SNRyconf2(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, SNRyconf3(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h3=plot(x, meanSNR(3,:), 'Color', [119 172 48]/255, 'LineStyle', '-', 'LineWidth', 2);            
end
ax = gca;
ax.Box = 'off';
ax.FontName = 'Times New Roman';
ax.GridLineStyle = '-.';
ax.GridColor = 'k';
ax.XGrid = 'on';
ax.YGrid = 'on';
ax.LineWidth = 1;
ax.XMinorTick = 'on';
ax.YMinorTick = 'on';
ax.TickDir = 'in';
ax.FontSize = fontsize-3;
ylabel('SNR', 'Interpreter','latex',FontSize=fontsize+2)
xlabel('$std(\mathbf{n})$', 'Interpreter','latex',FontSize=fontsize+2)
legend([h2, h1, h3],["DRUSvar", "DRUSmean", "$\mathbf{By}$"], ...
    "Interpreter","latex", "Orientation","vertical", "Location","best", ...
    fontsize=fontsize);
legend boxoff
xlim([0.018,max(x)]);


% save the charts as PDF
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
set(gcf, 'Renderer', 'opengl');
%print(gcf,[pwd '/numerical/Test_cysts/images/SNR'] ,'-dpdf','-r0')


