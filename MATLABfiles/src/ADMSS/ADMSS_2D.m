function Im_filt=ADMSS_2D(Im0,Mask,AD_param)
%ADMSS_2D - Two dimensional implementation of the Anisotropic Diffusion
%filter with Memory based on Speckle Statistics (ADMSS), recently
%published in 
%G. Ramos-Llorden et al., "Anisotropic Diffusion Filter With Memory Based 
%on Speckle Statistics for Ultrasound Images," IEEE Trans. Image Process., 
%vol.24, no.1, pp.345,358, Jan. 2015.
%
%-----DESCRIPTION-----
%ADMSS filters image 'Im0' by preserving speckle in relevant tissues and removing
%it in meaningless region. It makes use of a memory mechanism driven by
%Bayesian probabilty maps. Memory is switched off in meaningles region
%and it is activated in important tissues to preserve inital local
%structure tensor through the diffusion process.
%
%ADMSS implementation can work in two different modes: Bayesian probability maps of
%tissues are estimated a priori and provided or they are learned in situ in every
%iteration. Mode selection and the rest of ADMSS parameters are selected through the
%option structure 'AD_param' (See Inputs description).
%
%-----COPYRIGHT-----
%You can redistribute ADMSS_2D  and/or modify it under the terms of the GNU General 
%Public License as published by the Free Software Foundation, either version 
%3 of the License, or (at your option) any later version (attached in
%companion with this sofware) ADMSS_2D is provided on an "as is" basis, and
%the authors and distributors have no obligation to provide maintenance, 
%support, updates, enhancements, or modifications.
%
%If you use this sofware for research purpose, please cite
%
%G. Ramos-Llorden et al., "Anisotropic Diffusion Filter With Memory Based 
%on Speckle Statistics for Ultrasound Images," IEEE Trans. Image Process., 
%vol.24, no.1, pp.345,358, Jan. 2015.
%
%----HOW TO USE IT-----
% Syntax:  Im_filt=ADMSS_2D(Im0,AD_param)
%
% Inputs:
%    Im0 - [NxM] 2D Image to be filtered
%    Mask - [NxM] 2D Mask to avoid background (0s in background, 1s image)
%    AD_param - Structure with the following parameters:
%       AD_param.sigma  - Standard deviation sigma of Weickert's local
%       structure tensor Eq.(6)
%       AD_param.rho,  - Standard deviation rho of Weickert's local
%       structure tensor Eq.(6)
%       AD_param.nitmax - Maximum number of iterations
%       AD_param.n_memory - Polynomial degree in memory funct g() Eq.(21)
%       AD_param.delta_t - Time step
%       AD_param.estim - Option for the tissue probability estimation
%                        'Gaussian': Probabilities are estimated using a
%                        Gaussian Mixture Model
%                        'Gamma': Probabilities are estimated using a
%                        Gamma Mixture Model
%                        p_maps:  Already estimated probability maps [NxMxL]. 
%                               It is assummed that the last probability maps represents the meaningless tissue, blood for example,
%                        and hence memory should be switched off there
%    
%
% Outputs:
%    Im_filt - [NxM] 2D filtered image
%
% Example: 
%    Im0=double(imread([pwd,'/images/IVUS.png']));
%    AD_param.sigma=0.1;
%    AD_param.rho=0.1;
%    AD_param.nitmax=90;
%    AD_param.n_memory=15;
%    AD_param.delta_t=0.5; 
%    AD_param.estim='Gamma'; 
%    Im_filt=ADMSS_2D(Im0,AD_param);
%
% Other m-files required: DDE.m and those included in /Utils
%
%-----REFERENCES-----
%[1] G. Ramos-Llorden et al., "Anisotropic Diffusion Filter With Memory Based on Speckle Statistics for Ultrasound Images," IEEE Trans. Image Process., vol.24, no.1, pp.345,358, Jan. 2015
%[2] G. Vegas-Sanchez-Ferrero et al.," Probabilistic-driven oriented speckle reducing anisotropic diffusion with application to cardiac ultrasonic images.,", In Medical Image Computing and Computer-Assisted Intervention–MICCAI 2010 (pp. 518-525). Springer Berlin Heidelberg. 
%[3] G. Vegas-Sanchez-Ferrero et al., “Gamma mixture classifier for plaque detection in intravascular ultrasonic images,” IEEE Trans. Ultrason., Ferroelectr., Freq. Control, vol. 61, no. 1, pp. 44–61, Jan. 2014.
%[4] J. Weickert et al., “Efficient and reliable schemes for nonlinear diffusion filtering,” IEEE Trans. Image Process., vol. 7, no. 3, pp. 398–410, Mar. 1998.
%[5] M. A. Figueiredo and  A. K. Jain, “Unsupervised learning of finite mixture models,” IEEE Transactions on Pattern Analysis and Machine Intelligence, 24(3), 381-396.

%----AUTHORSHIP-----
% Code created by Gabriel Ramos Llordén and Gonzalo Vegas-Sánchez Ferrero
% Contact
% email: Gabriel.Ramos-LLorden@uantwerpen.be
% Website: http://visielab.uantwerpen.be/people/gabriel-ramos-llorden
% Jul 2015; Last revision: 27-Jul-2015


%% FLAGS and input detection
FLAG_ADMSS_WAIT=1;
FLAG_SELECTIVE_FILT=1;
FLAG_VOLTERRA=1;

if strcmp(AD_param.estim,'Gaussian')
    FLAG_EST=1;
    FLAG_GAUSS=1;
    FLAG_GAMMA=0;
elseif strcmp(AD_param.estim,'Gamma')
    FLAG_EST=1;
    FLAG_GAMMA=1;
    FLAG_GAUSS=0;
else
    FLAG_EST=0;
    p_maps=AD_param.estim;
    p_tissue=squeeze(p_maps(:,:,size(p_maps,3)));
end
if FLAG_ADMSS_WAIT,  hwait = waitbar(0,'ADMSS: Diffusing Image'); end
%% Parameters setting
delta_t=AD_param.delta_t;
sigma=AD_param.sigma;
rho=AD_param.rho;
n_memory=AD_param.n_memory;
nitmax=AD_param.nitmax;
dim=size(Im0);
im=Im0;
%% ADMSS starts to work here...
for nit = 1:nitmax
    if FLAG_EST  %If needed, parameter maps are re-estimated .
        r_estim=5; %Every r_estim iterations we re-estimate
        N_clas=3; %Number of classes in the Gamma mixture Model or Gaussian mixture model.
        if  mod(nit,r_estim)==0 || nit==1
            % PDF Mixture Model estimation Eq.(5) .
            if FLAG_GAMMA==1
              [prob, alpha_gamma, theta_gamma] = GMMestimator(im(Mask==1),N_clas,70,1e-5,0); %Speckle in IVUS is best modelled by Gammas [3] 
              mu=alpha_gamma.*theta_gamma;
              [~,I] = sort(mu,'descend'); %We re-order the parmeters so as to have in the last position, the tissue to be filtered (lowest mean, normally blood)
              p_maps = pmapsgamma(im, Mask, prob(I),alpha_gamma(I),theta_gamma(I)); %Gamma probability maps are created
            elseif FLAG_GAUSS==1
              [~,bestpp,bestmu,bestcov,~,~] = mixtures4(im(Mask==1)',N_clas,N_clas,0,1e-4,1); %Figueiredo's code for Gaussian [5];
              prob=bestpp';         
              mu=bestmu';
              vari=squeeze(bestcov)';
              [~,I] = sort(mu,'descend');
              p_maps = pmapsgauss(im,[],prob(I),mu(I),vari(I));  %Gaussian probability maps are created
            end
             p_tissue=p_maps(:,:,size(p_maps,3)); %Tissue to be removed Pg.350. In this case, we select the tissue with lowest mean (blood).
        end
    end 
    %% This block of code creates determines the diffusion tensor D.
    anterior = zeros(dim) - inf;
    v1x = zeros(dim);
    v1y=v1x; v2x=v1y; v2y=v2x;
    lambda_1 = ones(dim);
    lambda_2 = lambda_1;
    for i=1:size(p_maps,3)
        usigma=imgaussian(p_maps(:,:,i),sigma,4*sigma); %Gaussian filtered maps
        ux=derivatives(usigma,'x');
        uy=derivatives(usigma,'y');
        [Jxx, Jxy, Jyy] = Structuretensor(ux,uy,rho); %Structure tensor for the i-th p_map is defined Eq.(6)
        [~,mu2,U1,V1,U2,V2]=Eigenvectors(Jxx,Jxy,Jyy); %Eigenvalues and eigenvectors are derived.
        [valor, index] = max([anterior(:),mu2(:)],[],2); %
        anterior = valor;
        v1x(index==2) = U1(index==2);
        v1y(index==2) = V1(index==2);
        v2x(index==2) = U2(index==2);
        v2y(index==2) = V2(index==2); %We retain the eigenvectors where the i-th p_map prestend the maximum variation. Same philosophy as in [2]
        Grad=filtro_canny(p_maps(:,:,i),1*sigma); %edge filter for p_maps
        lambda_2(index==2) = 1-Grad(index==2); %If edge was high, esto==1, beta=0, anisotropic
        tolerdown=10^-5;
        tolerup=10^-1;
        lambda_2(lambda_2<tolerdown) = 0;
        lambda_2(lambda_2>(1 - tolerup)) = 1; %Tolerance thresolding to avoid unstabilities
    end
    % In case of no diagonalization, we start with an arbitrary basis
    index = sqrt(v1x.^2 + v1y.^2)<0.9;
    v1x(index) = 1;
    v1y(index) = 0;
    index = sqrt(v2x.^2 + v2y.^2)<0.9;
    v2x(index) = 0;
    v2y(index) = 1;
    index = (lambda_1==1 & lambda_2==1); %In anisotropic areas, i.e., lambda_1=lambda_2, we can increase the diffusivity deff
    deff=1.5*3;
    lambda_1(index==1)=deff;
    lambda_2(index==1)=deff;
    %% Memory mechanism: Volterra regularization. It determines the tensor L   
    if nit==1
        if FLAG_SELECTIVE_FILT==1
            lambda_1_filt=lambda_1.*p_tissue;
            lambda_2_filt=lambda_2.*p_tissue; %Selective filtering: Diffusion tensor S{D}. Eq.(20)
        else
            lambda_1_filt=lambda_1;
            lambda_2_filt=lambda_2;
        end
        Lxx = lambda_1_filt.*v1x.^2 + lambda_2_filt.*v2x.^2; %We recover the tensor L with the outer product formula Eq.(24)
        Lxy = lambda_1_filt.*v1x.*v1y + lambda_2_filt.*v2x.*v2y;
        Lyy = lambda_1_filt.*v1y.^2 + lambda_2_filt.*v2y.^2;
        gamma_1=lambda_1_filt;
        gamma_2=lambda_2_filt;
    else
        if FLAG_VOLTERRA==1
            [Lxx,Lxy,Lyy,gamma_1,gamma_2]=DDE(n_memory,lambda_1,lambda_2,v1x,v2x,v1y,v2y,p_tissue,gamma_1,gamma_2,FLAG_SELECTIVE_FILT); %Discretized Differential Delay equation Eq.(23)
        else
            Lxx = lambda_1.*v1x.^2 + lambda_2.*v2x.^2;
            Lxy = lambda_1.*v1x.*v1y + lambda_2.*v2x.*v2y;
            Lyy = lambda_1.*v1y.^2 + lambda_2_filt.*v2y.^2;
        end
    end
   %%
    im = diffusion_scheme_2D_implicit(im.*Mask,Lxx,Lxy,Lyy,delta_t); %2D Weickert semi implitic discretization [4]
    if FLAG_ADMSS_WAIT, waitbar(nit/nitmax,hwait); end
end
if FLAG_ADMSS_WAIT, close(hwait); end
Im_filt=im;
end

function [Lambda1,Lambda2,I2x,I2y,I1x,I1y]=Eigenvectors(Jxx,Jxy,Jyy)
% It gives the eigenvalues and eigenvector of the Structure Tensor whose
% components are Jxx, Jxy, Jyy (note that is symetric)
%
%  [Lambda1,Lambda2,I2x,I2y,I1x,I1y]=EigenVectors(Jxx,Jxy,Jyy)
%
% inputs, 
%   Jxx, Jxy and Jyy : Matrices with the values of the Hessian tensors
% 
% outputs,
%   Lambda1,Lambda2 : Eigen values
%   I2x,I2y,I1x,I1y : Eigen vectors
% 
% Compute the eigenvectors of J, v1 and v2
v2x = 2*Jxy; v2y = Jyy - Jxx + sqrt((Jxx - Jyy).^2 + 4*Jxy.^2);

% Normalize
normalization = sqrt(v2x.^2 + v2y.^2); 
i = (normalization ~= 0);
v2x(i) = v2x(i)./normalization(i);
v2y(i) = v2y(i)./normalization(i);

% Because eigenvectors are orthogonal 
v1x = -v2y; 
v1y = v2x;
% Compute the eigenvalues
mu1 = 0.5*(Jxx + Jyy + sqrt((Jxx - Jyy).^2 + 4*Jxy.^2));
mu2 = 0.5*(Jxx + Jyy - sqrt((Jxx - Jyy).^2 + 4*Jxy.^2));

% Sort eigen values by absolute value abs(Lambda1)<abs(Lambda2)
check=abs(mu1)>abs(mu2);
Lambda1=mu1; Lambda1(check)=mu2(check);
Lambda2=mu2; Lambda2(check)=mu1(check);
I1x=v1x; I1y=v1y; I2x=v2x; I2y=v2y; 
I1x(check)=v2x(check); I1y(check)=v2y(check);
I2x(check)=v1x(check); I2y(check)=v1y(check);
end

function [Jxx, Jxy, Jyy]=Structuretensor(ux,uy,rho)
% This function calculates the 2D
% regularized version of the structure tensor, i.e.,
% Jp ( grad(u) ) = conv ( Kp  , grad(u) * grad(u)^T ) 
% More information can be found at [4] or Weickert's PhD thesis.

% J(grad u_sigma)
Jxx = ux.^2;
Jxy = ux.*uy;
Jyy = uy.^2;
% Do the gaussian smoothing of the structure tensor
Jxx = imgaussian(Jxx,rho,6*rho);  %For small rho, convolution option is faster. If rho is big, we recommend implementing gaussian filters recursively
Jxy = imgaussian(Jxy,rho,6*rho);
Jyy = imgaussian(Jyy,rho,6*rho);
end