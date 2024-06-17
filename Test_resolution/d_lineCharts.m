%%% (scatterers) Display the line charts in Fig.3
%% evaluate and save the metrics
clc
close all
clear

flag_display = 0;
flag_avg = 0;  % taking the mean of the scores of ROIS in one image
channels = 1;
addpath(genpath('src'));
path_phantom = 'phantom_9.hdf5';
imgSize = 256;
scan = linear_scan(linspace(-0.018,0.018,imgSize).', linspace(0.01,0.036+0.01,imgSize).');
c = 11;  % number of levels of the additive noise
repeat = 20;
REPEAT = 20;
folder = '0.0-0.1/sigma1.2/'; 
m = [4,5,6,7,8,9,10,11,12]; % 9 realizations of the multiplicative noise 
filepath = [pwd '/Test_resolution/' folder];
metrics = cell(length(m),4,7);


for p = 1:length(m)
    filename = ['sigma12_m' num2str(m(p))];
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
    image.postenv = 1;
    [yFWHMA(i,:), yFWHML(i,:), yCNR(i), ygCNR(i), ySNR(i), yKS(i)] = evaluation(path_phantom, image, flag_display, flag_avg); 
    image.postenv = 2;
    [~, yFWHML(i,:), yCNR(i), ygCNR(i), ySNR(i), yKS(i)] = evaluation(path_phantom, image, flag_display, flag_avg); 
    end
    
    % Mean
    for i = 1: c
    image = us_image(); image.scan = scan; image.number_plane_waves=1;
    tmp = abs(squeeze(meanImgs(i, :, :)));
    image.data = tmp./max(tmp(:)); 
    image.postenv = 1;
    [MeanFWHMA(i,:), MeanFWHML(i,:), MeanCNR(i), MeangCNR(i), MeanSNR(i), MeanKS(i)] = evaluation(path_phantom, image, flag_display ,flag_avg); 
    image.postenv = 2;
    [~, MeanFWHML(i,:), MeanCNR(i), MeangCNR(i), MeanSNR(i), MeanKS(i)] = evaluation(path_phantom, image, flag_display ,flag_avg); 
    end
    
    % var
    for i = 1: c
    image = us_image(); image.scan = scan; image.number_plane_waves=1;
    tmp = reshape(abs(squeeze(varImgs(i, :, :))),[imgSize,imgSize]);
    image.data = tmp./max(tmp(:)); 
    image.postenv = 1;
    [VarFWHMA(i,:), VarFWHML(i,:), VarCNR(i), VargCNR(i), VarSNR(i), VarKS(i)] = evaluation(path_phantom, image, flag_display ,flag_avg);
    image.postenv = 2;
    [~, VarFWHML(i,:), VarCNR(i), VargCNR(i), VarSNR(i), VarKS(i)] = evaluation(path_phantom, image, flag_display ,flag_avg);
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
% save([pwd '/Test_resolution/0.0-0.1/metrics.mat'], 'metrics' 'c'); 

%% load the metrics and plot
clc
load("Test_resolution/0.0-0.1/metrics.mat")
type = 1;  fontsize=19;

[ml, ~, ~] = size(metrics);
% 25 scatterers in one multiplicative noise realization
% ml=9 realizations of multiplicative noise
% 3 types of images (the degraded one , DRUSmean, DRUSvar)
% c = 11 levels of additive noise
A = zeros(25*ml,3,c);  
L = zeros(25*ml,3,c);
for p = 1:ml
    for q = 1: 25
    A(25*(p-1)+q, 1,:) = metrics{p,1,1}(:,q);
    A(25*(p-1)+q, 2,:) = metrics{p,2,1}(:,q);
    A(25*(p-1)+q, 3,:) = metrics{p,4,1}(:,q);
    L(25*(p-1)+q, 1,:) = metrics{p,1,2}(:,q);
    L(25*(p-1)+q, 2,:) = metrics{p,2,2}(:,q);
    L(25*(p-1)+q, 3,:) = metrics{p,4,2}(:,q);
    end
end
meanA = squeeze(mean(A)); meanL = squeeze(mean(L));
stdA = squeeze(std(A)); stdL = squeeze(std(L));

% compute the bounds of the std
AlowerBound1 = flip((meanA - stdA),2);
AupperBound1 = meanA + stdA;
AlowerBound2 = flip((meanA - 2 * stdA),2);
AupperBound2 = meanA + 2 * stdA;
AlowerBound3 = flip((meanA - 3 * stdA),2);
AupperBound3 = meanA + 3 * stdA;

LlowerBound1 = flip((meanL - stdL),2);
LupperBound1 = meanL + stdL;
LlowerBound2 = flip((meanL - 2 * stdL),2);
LupperBound2 = meanL + 2 * stdL;
LlowerBound3 = flip((meanL - 3 * stdL),2);
LupperBound3 = meanL + 3 * stdL;

% prepare the x-axis coordinates for drawing
x =  [0,0.003,0.006,0.009,0.012,0.015,0.018,0.04,0.06,0.08,0.1];
xconf = [x flipud(x')'];
Ayconf1 = [AupperBound1 AlowerBound1];
Ayconf2 = [AupperBound2 AlowerBound2];
Ayconf3 = [AupperBound3 AlowerBound3];

Lyconf1 = [LupperBound1 LlowerBound1];
Lyconf2 = [LupperBound2 LlowerBound2];
Lyconf3 = [LupperBound3 LlowerBound3];

% chart of FWHM_A
figure('Position', [476 356 700 370]);hold on;
switch type    
    case 1        
        fill(xconf, Ayconf1(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h1=plot(x, meanA(1,:), 'Color', [0 72 186]/255, 'LineStyle', '-', 'LineWidth', 2);  
        plot(x,meanA(1,:),'o','MarkerSize',8,'MarkerFaceColor',[0 72 186]/255,'MarkerEdgeColor',[0 72 186]/255,'HandleVisibility', 'off')        

        fill(xconf, Ayconf1(2,:), [219, 112, 147] / 255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h2=plot(x, meanA(2,:), 'Color', [219, 112, 147] / 255, 'LineStyle', '-', 'LineWidth', 2); 
        plot(x,meanA(2,:),'o','MarkerSize',8,'MarkerFaceColor',[219, 112, 147] / 255,'MarkerEdgeColor',[219, 112, 147] / 255,'HandleVisibility', 'off')        

        fill(xconf, Ayconf1(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h3=plot(x, meanA(3,:), 'Color', [119 172 48]/255, 'LineStyle', '-', 'LineWidth', 2);  
        plot(x,meanA(3,:),'o','MarkerSize',8,'MarkerFaceColor',[119 172 48]/255,'MarkerEdgeColor',[119 172 48]/255,'HandleVisibility', 'off')        

    case 2        
        fHdl(1) = fill(xconf, Ayconf1(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, Ayconf2(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, Ayconf3(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h1=plot(x, meanA(1,:), 'Color', [0 72 186]/255, 'LineStyle', '-', 'LineWidth', 2);    
        fHdl(1) = fill(xconf, Ayconf1(2,:), [219, 112, 147] / 255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, Ayconf2(2,:), [219, 112, 147] / 255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, Ayconf3(2,:), [219, 112, 147] / 255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h2=plot(x, meanA(2,:), 'Color', [219, 112, 147] / 255, 'LineStyle', '-', 'LineWidth', 2);  
        fHdl(1) = fill(xconf, Ayconf1(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, Ayconf2(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, Ayconf3(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h3=plot(x, meanA(3,:), 'Color', [119 172 48]/255, 'LineStyle', '-', 'LineWidth', 2);            
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
ylabel('FWHM(axial) [mm]', 'Interpreter','latex',FontSize=fontsize+2)
xlabel('$std(\mathbf{n})$', 'Interpreter','latex',FontSize=fontsize+2)
legend([h2, h1, h3],["DRUSvar", "DRUSmean", "$\mathbf{By}$"], ...
    "Interpreter","latex", "Orientation","vertical", "Location","best", ...
    fontsize=fontsize);
legend boxoff
ylim([0.38,1.25]);
xlim([0.018,max(x)]);


% chart of FWHM_L
figure('Position', [476 356 700 370]);hold on;
switch type    
    case 1        
        fill(xconf, Lyconf1(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h1=plot(x, meanL(1,:), 'Color', [0 72 186]/255, 'LineStyle', '-', 'LineWidth', 2, 'DisplayName','Mean');  
        plot(x,meanL(1,:),'o','MarkerSize',8,'MarkerFaceColor',[0 72 186]/255,'MarkerEdgeColor',[0 72 186]/255,'HandleVisibility', 'off')        

        fill(xconf, Lyconf1(2,:), [219, 112, 147] / 255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h2=plot(x, meanL(2,:), 'Color', [219, 112, 147] / 255, 'LineStyle', '-', 'LineWidth', 2,'DisplayName','Var'); 
        plot(x,meanL(2,:),'o','MarkerSize',8,'MarkerFaceColor',[219, 112, 147] / 255,'MarkerEdgeColor',[219, 112, 147] / 255,'HandleVisibility', 'off')        

        fill(xconf, Lyconf1(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h3=plot(x, meanL(3,:), 'Color', [119 172 48]/255, 'LineStyle', '-', 'LineWidth', 2,'DisplayName','y');  
        plot(x,meanL(3,:),'o','MarkerSize',8,'MarkerFaceColor',[119 172 48]/255,'MarkerEdgeColor',[119 172 48]/255,'HandleVisibility', 'off')        

    case 2        
        fHdl(1) = fill(xconf, Lyconf1(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, Lyconf2(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, Lyconf3(1,:), [0 72 186]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h1=plot(x, meanL(1,:), 'Color', [0 72 186]/255, 'LineStyle', '-', 'LineWidth', 2);    
        fHdl(1) = fill(xconf, Lyconf1(2,:), [219, 112, 147] / 255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, Lyconf2(2,:), [219, 112, 147] / 255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, Lyconf3(2,:), [219, 112, 147] / 255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h2=plot(x, meanL(2,:), 'Color', [219, 112, 147] / 255, 'LineStyle', '-', 'LineWidth', 2);  
        fHdl(1) = fill(xconf, Lyconf1(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        fHdl(2) = fill(xconf, Lyconf2(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);       
        fHdl(3) = fill(xconf, Lyconf3(3,:), [119 172 48]/255, 'EdgeColor', 'none', 'FaceAlpha', .2);        
        h3=plot(x, meanL(3,:), 'Color', [119 172 48]/255, 'LineStyle', '-', 'LineWidth', 2);            
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
ylabel('FWHM(lateral) [mm]', 'Interpreter','latex',FontSize=fontsize+2)
xlabel('$std(\mathbf{n})$', 'Interpreter','latex',FontSize=fontsize+2)
legend([h2, h1, h3],["DRUSvar", "DRUSmean", "$\mathbf{By}$"], ...
    "Interpreter","latex", "Orientation","vertical", "Location","best", ...
    fontsize=fontsize);
legend boxoff
xlim([0.018,max(x)]);

%% save as PDF
set(gcf,'Units','Inches');
pos = get(gcf,'Position');
set(gcf,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(gcf,['Test_resolution/images/-6dB_Lateral_y'] ,'-dpdf','-r0')



