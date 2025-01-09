function Maps = pmapsgamma(I0, Mask, W, K,THETA)
%this function receives the parameters of a Gamma Mixture Model and
%creates the probability maps%
origen = 1;
destino = length(W);
denominador = zeros(size(I0));
dim = size(I0);
for i=1:length(W)
    denominador = denominador + W(i)*gampdf(I0,K(i),THETA(i)); %fx
end
Maps = zeros([size(I0),length(origen:destino)]);
for i=origen:length(W)
    numerador = (W(i)*gampdf(I0,K(i),THETA(i)));
    esto = numerador./(denominador+eps);
    esto(isnan(esto)) = 1;
    Maps(:,:,i) = esto;
end
Maps(Maps<0) = 0;
Maps(Maps>1) = 1;
