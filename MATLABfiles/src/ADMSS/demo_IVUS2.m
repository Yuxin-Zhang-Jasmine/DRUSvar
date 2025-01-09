%In this demo, an IVUS image is filtered with the Anisotropic Diffusion
%Filter with memory. Probabilty maps are learned internally in each
%iteration.
%
% Original image was obtained 
%from https://www.osapublishing.org/oe/fulltext.cfm?uri=oe-16-16-12313&id=170288
%%
clear all, clc, close all
addpath(strcat(pwd,'/utils'));
Im0=double(imread([pwd,'/images/IVUS.png']));
%% Parameters definition
AD_param.sigma=0.1;
AD_param.rho=0.1;
AD_param.nitmax=30;
AD_param.n_memory=15;
AD_param.delta_t=0.5;
AD_param.estim='Gamma';
%% Filtering
Im_filt=ADMSS_2D(Im0,ones(size(Im0)),AD_param); %No background. Mask contains only 1s
%% Visualization
disp('Done!, now visualizing results');
figure(1)
imshow([Im0,Im_filt],[])
title('Original                          Filtered')
