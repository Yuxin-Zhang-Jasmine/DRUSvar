% GMMestimator - implementation of the Gamma Mixture Model
% 
% -----description-----
% GMMestimator estimates the components of a finite
% mixture model following a Gamma distribution.
% 
% -----copyright-----
% You can redistribute GMMestimator and/or modify it under the terms of the gnu general 
% public license as published by the free software foundation, either version 
% 3 of the license, or (at your option) any later version (attached in
% companion with this sofware) GMMestimator is provided on an "as is" basis, and
% the authors and distributors have no obligation to provide maintenance, 
% support, updates, enhancements, or modifications.
% 
% If you use this sofware for research purpose, please cite
% 
% G. Vegas-Sánchez-Ferrero, M. Martín-Fernández, J. Miguel Sanches "A Gamma
% Mixture Model for IVUS Imaging", Multi-Modality Atherosclerosis Imaging
% and Diagnosis. Editors: Luca Saba, João Miguel Sanches, Luís Mendes
% Pedro, Jasjit S. Suri ISBN: 978-1-4614-7424-1 (Print) 978-1-4614-7425-8
% (Online), Pages 155-171. 2014.
% 
% ----how to use it-----
% Syntax:  [w, alpha, beta] = GMMestimator(y,nl,maxIter,tol_error,flag_pinta,w_0,alpha_0,beta_0)
% Inputs:
%    y - vector of samples
%    nl - number of mixture components.
%    maxIter - Maximum number of iterations
%    tol_error - Tolerance assumed for convergence
%    flag_pinta - Flag to show the evolution of fitting
%    w_0 = initial weights size (1 x nl). They should sum up to 1. (optional)
%    alpha_0 = initial alpha parameter of each Gamma component (size: 1 x nl) (optional)
%    beta_0  = initial beta parameter of each Gamma component (size: 1 x nl) (optional)
%    
% outputs:
%    w - vector of weights for component
%    alpha - vector of alphas
%    beta - vector of betas
%    The components are ordered by decresing weights.
% 
% example: 
%    Y = [gamrnd(2,3,1,1000) gamrnd(2,5,1,1000) gamrnd(2,7,1,3000)];
%    [prob, alpha_gamma, beta_gamma] = GMMestimator(Y,3,70,1e-5,1);
% 
% other m-files required: dde.m and those included in /utils
% 
% -----references-----
% G. Vegas Sánchez-Ferrero, A. Tristán Vega, S. Aja-Fernández, M.
% Martín-Fernández, C. Palencia, R. Deriche. “Anisotropic LMMSE denoising
% of MRI based on statistical tissue models”. Symposium on Biomedical
% Imaging: From nano to macro (ISBI). Barcelona (Spain). 2012.
%
% Gonzalo Vegas-Sánchez-Ferrero, Santiago Aja Fernández, Marcos Martín
% Fernández, Alejandro Frangi, César Palencia. “Probabilistic-Driven
% Oriented Speckle Reducing Anisotropic Diffusion with Application to
% Cardiac Ultrasonic Images”. International Conference on Medical Image
% Computing and Computer Assisted Intervention (MICCAI), Beiging, Sep 2010.
%
% Gonzalo Vegas Sánchez-Ferrero, José  Seabra, Santiago Aja-Fernández,
% Marcos Martín-Fernández, César Palencia, Joao Sanches. "Gamma Mixture
% Classifier for Plaque Detection in Intravascular Ultrasonic Images". IEEE
% Transactions on Ultrasonics, Ferroelectrics, and Frequency Control.
% vol.61, no.1, pp.44,61, January 2014.
%
% Gonzalo Vegas-Sánchez-Ferrero, Santiago Aja-Fernández, César Palencia,
% Marcos Martín-Fernández. "Generalized Gamma Mixture Model for Ultrasonic
% Tissue Characterization". Computational and Mathematical Methods in
% Medicine. Article ID 481923, 25. 2013.
%
% H. Curiale, A. Haak, G. Vegas-Sánchez-Ferrero, B. Ren, S. Aja-Fernández,
% J. G. Bosch, "Fully automatic detection of salient features in 3D
% transesophageal images", Ultrasound in Medicine and Biology.
% 40(12):2868-84, 2014.
% 
% A. Haak, G. Vegas-Sánchez-Ferrero, H. W. Mulder, B. Ren, H. A. Kirisli,
% C. Metz, G. van Burken, M. van Stralen, J. P. W: Pluim, A. F.W. van der
% Steen, T. van Walsum, J. G. Bosch, "Segmentation of multiple heart
% cavities in 3D transesophageal ultrasound images", IEEE Transactions on
% Ultrasonics, Ferroelectrics, and Frequency Control. Vol. 62, No. 6, June
% 2015.
% 
% G. Ramos-Llordén, G. Vegas-Sánchez-Ferrero, M. Martín-Fernández, S.
% Aja-Fernández, C. Alberola-López, "Anisotropic diffusion filter with
% memory based on speckle statistics for medical images". IEEE Transactions
% on Image Processing. Vol. 24, No. 1, January 2015.
%
% A. Haak, B. Ren, H. W. ulder, G. Vegas-Sánchez-Ferrero, G. van Bruken, A.
% F.W. van der Steen, M. van Stralen, P.W. Pluim, T. van Walsum, J.G.
% Bosch,  “Improved segmentation of multiple cavities of the heart in
% wide-view 3D transesophageal echocardiograms”, Ultrasound in Medicine and
% Biology. Vol. 41, Issue 7 , 1991-2000, March 2015.
%
%
% ----authorship-----
% code created by Gonzalo Vegas-Sanchez-Ferrero
% contact
% email: gvegsan@lpi.tel.uva.es
% website: https://www.lpi.tel.uva.es/gvegsan
% sep 2015; last revision: 09-sep-2015

function [pesos, alpha, beta] = GMMestimator(y,nl,maxIter,tol_error,flag_pinta,pesos_inicial,alpha_inicial,beta_inicial)

y = y(:)';
scalex = sum(y)/length(y);
% scalex = 1;
y = y./scalex;

% flag_pinta = 0;

[N,X] = hist(y,linspace(min(y(:)),max(y(:)),300));

h = X(2)-X(1);
Y = N/(sum(N*h));

if nargin<6
    nl = min(floor(length(y)/10),nl);
    
    % Deterministic initialization.
    %     intensidad_min = min(y);
    %     intensidad_max = max(y);
    %     index = ones(size(y));
    %     paso = (intensidad_max-intensidad_min)/nl;
    %     for conta=1:nl-1
    %         index = index + (y>=paso*conta);
    %     end
    
    % Initialization with k-means. When void clusters, one can use the
    % previous deterministic initialization    
    index = kmeans(y',nl);
    
    pesos = zeros(1,nl); alpha = zeros(1,nl); beta = zeros(1,nl);
    for conta=1:nl
        pesos(conta) = sum(index==conta)/sum(N);        
        beta(conta) = var(y(index==conta))/mean(y(index==conta));
        alpha(conta) = mean(y(index==conta))/beta(conta);
    end
    
else
    pesos = pesos_inicial;
    alpha = alpha_inicial;
    beta = beta_inicial;
end

if flag_pinta
    figure(2);
    plot(X,Y,'k--');
    hold on;
    total2 = zeros(size(X));
    for conta = 1:nl
        plot(X,pesos(conta)*gampdf(X,alpha(conta),beta(conta)),'m');
        total2 = total2 + pesos(conta)*gampdf(X,alpha(conta),beta(conta));
    end
    plot(X,total2,'m','linewidth',2);
end


p = zeros(nl,length(y));
K = zeros(1,length(y));
for conta = 1:nl        
    p(conta,:) = gampdf(y,alpha(conta),beta(conta));    
    K = K + p(conta,:).*pesos(conta);     
end

% Weights calculation
w_t = zeros(1,length(y));
for n=1:nl,
    w_t(n,:) = p(n,:) .* pesos(n) ./ K;
end
% In case inf/inf, we set it as 1.
w_t(isnan(w_t)) = 1;

% Expectation Maximization (EM) algorithm: 
err = Inf;
i=1;
iter = i;
options = [];
options = statset(statset('gamfit'), options);

pesos0 = pesos;
alpha0 = alpha;
beta0 = beta;

lkeqn = @(x_lkeqn,y_lkeqn) y_lkeqn - log(x_lkeqn+eps) + psi(x_lkeqn+eps);

while (err>tol_error && i<maxIter),
    K = zeros(1,length(y));
    
    for n=1:nl,
        % Cálculo de los pesos
        pesos(n) = sum(w_t(n,:))/(sum(w_t(:))); 
                
        % Cáculo del parámetro alpha        
        
            % Definimos la constante        
            const1 = log(sum(w_t(n,:).*y)/sum(w_t(n,:)));
            const2 = -sum(w_t(n,:).*log(y+eps))/(sum(w_t(n,:))+eps);        
            const = const1 + const2;

            if lkeqn(alpha(n), const) > 0
                upper = alpha(n); lower = .5 * upper;
                while lkeqn(lower, const) > 0
                    upper = lower;
                    lower = .5 * upper;
                    if lower < realmin(class(y)) % underflow, no positive root
                        error('stats:gamfit:NoSolution',...
                            'Unable to reach a maximum likelihood solution.');
                    end
                end
            else
                lower = alpha(n); upper = 2 * lower;
                while lkeqn(upper, const) < 0
                    lower = upper;
                    upper = 2 * lower;
                    if upper > realmax(class(y)) % overflow, no finite root
                        error('stats:gamfit:NoSolution',...
                            'Unable to reach a maximum likelihood solution.');
                    end
                end
            end
            bnds = [lower upper];        

          try
            [ahat, lkeqnval, err] = fzero(lkeqn, bnds, options, const);
            alpha(n) = ahat;
          catch
%              keyboard;            
          end
              
        % Beta Calculation
        beta(n) = sum(w_t(n,:).*y)/(sum(w_t(n,:))*alpha(n)+eps);
        
        if isnan(beta(n)), keyboard, end
        if isnan(alpha(n)), keyboard, end
                
        p(n,:) = gampdf(y,alpha(n),beta(n));
        K = K + p(n,:).*pesos(n);         
    end               
        
            
    % Updates:
    w_t = zeros(1,length(y));       
    
    for n=1:nl,
        w_t(n,:) = p(n,:) .* pesos(n) ./ K;
        w_t(isnan(w_t)) = 1;
    end    
    
    err = 0;
    for n=1:nl, % check whether estimates converged or are still varying
        err = err + max(abs(pesos - sum(w_t,2)'/(sum(w_t(:)))));
    end
    if isnan(err), keyboard, end
    i = i+1;    
    iter = [iter, i];       
            
    if flag_pinta
        figure(3);
        cla
        plot(X,Y,'k--');
        hold on;
        total2 = zeros(size(X));
        total = zeros(size(X));
        for conta = 1:nl
            plot(X,pesos(conta)*gampdf(X,alpha(conta),beta(conta)),'r');
            plot(X,pesos0(conta)*gampdf(X,alpha0(conta),beta0(conta)),'b');
            
            total = total + pesos0(conta)*gampdf(X,alpha0(conta),beta0(conta));
            total2 = total2 + pesos(conta)*gampdf(X,alpha(conta),beta(conta));
        end
        plot(X,total,'b','linewidth',2);
        plot(X,total2,'r','linewidth',2);
        pause
    end
end

for n=1:nl
    pesos(n) = sum(w_t(n,:))/(sum(w_t(:)));
end
beta = beta*scalex;

[dum,idx1] = sort(pesos,'descend');
pesos = dum;            % coefficients
beta = beta(idx1);      % Gamma parameters
alpha = alpha(idx1);    % Gamma parameters

