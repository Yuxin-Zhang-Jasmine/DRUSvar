%%% evaluate the PICMUS in-vitro results (TABLE I)
close all
clear
clc

% Parameter(s)
model = 'DRUS'; % 'Deno'

% Other configuration
addpath(genpath('src'));
MulPath = '../3_repeatNumSensitivity/50it/';
flag_display = 1;
flag_avg = 1;
numPhantoms  = 4;
repeatIdx = 8;  % corresponds to repeat 10 times (50it x 10)
mapnum = 4;     % "median"; "var"; "mean"; "std"
if strcmp(model, 'Deno')
    mapIdxs = [3];       % [mean]
    filename = [model, 'MetricsmapEnv.mat'];
    load([MulPath filename])
elseif strcmp(model, 'DRUS')
    mapIdxs = [3, 2];     % [mean, var]
    filename = [model, 'MetricsmapEnv.mat'];
    load([MulPath filename])
elseif strcmp(model, 'ADMSS')
    mapIdxs = [3];       % [mean]
    filename = 'DRUSMetricsmapEnv.mat';
    load([MulPath filename])
    AD_param.sigma=4;
    AD_param.rho=4;
    AD_param.nitmax=31;
    AD_param.n_memory=9;
    AD_param.delta_t=0.3;
    AD_param.estim='Gamma';
    scan = linear_scan(linspace(-0.018,0.018,256).', linspace(0.01,0.036+0.01,256).');
    mask = ones(scan.Nz, scan.Nx);
end
phantomIdxs = [3,4]; % [ER, EC]

% Metrics matrix Initialization
FWHMAsmap = zeros(mapnum, numPhantoms);
FWHMLsmap = zeros(mapnum, numPhantoms);
CNRsmap  = zeros(mapnum, numPhantoms);
gCNRsmap = zeros(mapnum, numPhantoms);
SNRsmap  = zeros(mapnum, numPhantoms);
KSsmap   = zeros(mapnum, numPhantoms);

% Evaluation
for phantomIdx = phantomIdxs
    pathPhantom = ['picmus_phantom_' num2str(phantomIdx) '.hdf5'];
    for mapIdx = mapIdxs  
        Im = ImageSetsmap{repeatIdx}{mapIdx, phantomIdx};
        if strcmp(model, 'ADMSS') %-- ADMSS filtering
            Im.data = abs(ADMSS_2D(Im.data,mask,AD_param));
            Im.postenv=[];
        end
        [FWHMAsmap(mapIdx, phantomIdx), FWHMLsmap(mapIdx, phantomIdx), CNRsmap(mapIdx, phantomIdx), gCNRsmap(mapIdx, phantomIdx), SNRsmap(mapIdx, phantomIdx), KSsmap(mapIdx, phantomIdx)] = evaluation(pathPhantom, Im, flag_display, flag_avg);     
    end
end



%% [mean, var]
mapNum = length(mapIdxs); 

ERFWHMA = zeros(mapNum,1);
ERFWHML = zeros(mapNum,1);
ERgCNR  = zeros(mapNum,1);

ECSNR   = zeros(mapNum,1);
ECgCNR  = zeros(mapNum,1);

for i = 1:mapNum
    mapIdx = mapIdxs(i);
    % ER
    ERFWHMA(i) = FWHMAsmap(mapIdx, 3);
    ERFWHML(i) = FWHMLsmap(mapIdx, 3);
    ERgCNR(i) = gCNRsmap(mapIdx, 3);
    % EC
    ECgCNR(i) = gCNRsmap(mapIdx, 4);
    ECSNR(i) = SNRsmap(mapIdx, 4);
end
clearvars -except ERFWHMA ERFWHML ERgCNR ECgCNR ECSNR 