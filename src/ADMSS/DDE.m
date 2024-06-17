function [Lxx,Lxy,Lyy,gamma_1,gamma_2]=DDE(n_memory,lambda_1,lambda_2,v1x,v2x,v1y,v2y,p_tissue,gamma_1,gamma_2,FLAG_SELECTIVE_FILT)
%This function implements the discretization of the Volterra equation.
    Beta = (1- p_tissue)./(p_tissue.^n_memory + eps);  %Memory function Eq.[21]
    if FLAG_SELECTIVE_FILT==1
        lambda_1_filt=lambda_1.*p_tissue;
        lambda_2_filt=lambda_2.*p_tissue; %Selective filtering: Diffusion tensor S{D}.
    else
        lambda_1_filt=lambda_1;
        lambda_2_filt=lambda_2;
    end
    gamma_1 = 1./(1+Beta).*(Beta.*gamma_1  + lambda_1_filt);
    gamma_2 = 1./(1+Beta).*(Beta.*gamma_2  +  lambda_2_filt);
    Lxx = gamma_1.*v1x.^2   + gamma_2.*v2x.^2;
    Lxy = gamma_1.*v1x.*v1y   + gamma_2.*v2x.*v2y;
    Lyy = gamma_1.*v1y.^2   + gamma_2.*v2y.^2;