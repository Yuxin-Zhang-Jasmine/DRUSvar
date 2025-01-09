%In this demo, a 2D echocardiography sequence is filtered with the Anisotropic Diffusion
%Filter with memory. Probabilty maps are learned a priori
%%
clear all, clc, close all
addpath(strcat(pwd,'/utils'));
load([pwd,'/images/Seq.mat']);
Mask_big=ones(size(Seq));
Mask=~(squeeze(Seq(:,:,1))==0);
Mask_big=bsxfun(@times,Mask,Mask_big);
%% A priori probability map estimation. This helps to make the filtering process faster...
N_classes=3;
FLAG_STAY=1;
if exist('p_maps')==0
       while (FLAG_STAY)
            [bestk,bestpp,bestmu,bestcov,dl,countf] = mixtures4(Seq(Mask_big==1)',N_classes,N_classes,0.01,1e-6,1); %Gaussian Mixture Model has been assumed.
             Error=dl(end);
             if Error<(1e+6)*3.78, FLAG_STAY=0; end
       end
       disp(sprintf('Probability maps have been learned. Now filtering...\n'))
       prob=bestpp';
       mu=bestmu';
       vari=squeeze(bestcov);
       [B,indx] = sort(mu,'descend');
end
%% Parameter definition
AD_param.sigma=0.1;
AD_param.rho=0.1;
AD_param.nitmax=40;
AD_param.n_memory=17;
AD_param.delta_t=0.5;
%% Filtering
N_tmax=size(Seq,3);
for N_t=1:N_tmax
    I0=squeeze(Seq(:,:,N_t));
    p_maps = pmapsgauss(I0,[],prob(indx),mu(indx),vari(indx));
    AD_param.estim=p_maps;
    Im_filt=ADMSS_2D(I0,Mask,AD_param);
    Seq_filt(:,:,N_t)=Im_filt;
    pause(2)
    disp(sprintf('The slice number %d has been processed.', N_t))
end
disp('Done!, now visualizing results');
%% Visualization
FLAG_view=1;
N_t=1;
Minint=min(Seq(:));
Maxint=max(Seq(:));
while(FLAG_view)
   figure(1)
   imshow([squeeze(Seq(:,:,N_t)),squeeze(Seq_filt(:,:,N_t))],[Minint,Maxint]);
   title('Original                          Filtered')
   drawnow
   N_t=N_t+1;
   pause(1/200)
   if N_t==(N_tmax+1) ,N_t=1; end
end
