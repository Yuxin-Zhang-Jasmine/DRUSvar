function Maps = pmapsgauss(I0,Mask,prob,mu,vari)
%this function receives the parameters of a Gaussian Mixture Model and
%creates the probability maps%
N_maps=length(prob);
dim = size(I0);
I=I0;
fx=zeros(dim);
for i=1:N_maps
    fx = fx + prob(i)*normpdf(I,mu(i),sqrt(vari(i)));
end
Maps = zeros([size(I0),N_maps]);
for i=1:N_maps
    numerador = prob(i)*normpdf(I,mu(i),sqrt(vari(i)));
    esto = numerador./(fx+eps);
    esto(isnan(esto)) = 1;
    Maps(:,:,i) = esto;
end
Maps(Maps<0) = 0;
Maps(Maps>1) = 1; %probability maps are always bounded [0,1]
