%In this demo, a 2D Ultrasound image is filtered with the Anisotropic Diffusion
%Filter with memory. Probabilty maps are learned in each iteration
%
%G. Ramos-Llorden et al., "Anisotropic Diffusion Filter With Memory Based 
%on Speckle Statistics for Ultrasound Images," IEEE Trans. Image Process., 
%vol.24, no.1, pp.345,358, Jan. 2015.

%Original image obtained from 
%http://telin.ugent.be/~sanja/Sanja_files/UltrasoundDemo.htm
%%
clear, clc, close all
addpath(strcat(pwd,'/utils'));
Im0=double(rgb2gray(imread([pwd,'/images/US.jpg'])));
Mask=double((imread([pwd,'/images/Mask.png'])));
%% Parameters definition
AD_param.sigma=0.8;
AD_param.rho=0.8;
AD_param.nitmax=55;
AD_param.n_memory=25;
AD_param.delta_t=0.1;
AD_param.estim='Gaussian';
%% Filtering
Im_filt=ADMSS_2D(Im0,Mask,AD_param);
%% Visualization
disp('Done!, now visualizing results');
figure(1)
imagesc([Im0.*Mask,Im_filt.*Mask])
axis off
colormap gray
title('Original                          Filtered')