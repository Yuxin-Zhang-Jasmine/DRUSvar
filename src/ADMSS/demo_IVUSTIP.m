%In this demo, a IVUS image in polar form is filtered with the Anisotropic Diffusion
%Filter with memory. 
%
%%
clear all, clc, close all
addpath(strcat(pwd,'/utils'));
load(strcat(pwd,'/images/rfEnvmaps.mat')) %In this case, p_maps have been learned a priori
load(strcat(pwd,'/images/rfEnv.mat')) %The image comes in polar form, (rows-->radio, columns--> angle). 
Im0=rfEnv;
%% Parameters definition
AD_param.sigma=0.05;
AD_param.rho=0.05;
AD_param.nitmax=80;
AD_param.n_memory=5;
AD_param.delta_t=0.15;
AD_param.estim=p_maps;
%% Filtering
Im_filt=ADMSS_2D(Im0,ones(size(Im0)),AD_param); %No background, mask contains only 1s
%% Visualization
disp('Done!, now visualizing results');
Im_filt=polar2cart(Im_filt,256,inf,10); %We transform the polar image into Cartesian form [256x256]
Im0=polar2cart(Im0,256,inf,10);
figure(1)
imshow(log(1+abs([Im0,Im_filt])),[]) %We apply compression for better visualization
title('Original                          Filtered')