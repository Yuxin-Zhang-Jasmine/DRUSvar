% Note: degradation is applied directly in python 
% the plot here is used to have an intuition

sigma = 1.2;  % corresponding to 0.17 mm
kerHalfLen = ceil(2 * sigma);
x = -kerHalfLen:1:kerHalfLen;
ycos = cos(2*pi*0.5*x);
ygauss = exp(-0.5 .* (x./sigma).^2);

y = ycos .* ygauss;
figure;
subplot(3,1,1); plot(x, ycos)
subplot(3,1,2); plot(x, ygauss);
subplot(3,1,3); plot(x, y);