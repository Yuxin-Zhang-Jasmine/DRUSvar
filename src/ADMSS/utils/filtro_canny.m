function mag=filtro_canny(a,sigma);

m = size(a,1);
n = size(a,2);
rr = 2:m-1; cc=2:n-1;

e = repmat(false, m, n);

Sigma = 1.0;          % Default Std dev of gaussian for canny
  threshSpecified = 0;  % Threshold is not yet specified

% Magic numbers
  GaussianDieOff = .0001;  
  PercentOfPixelsNotEdges = .7; % Used for selecting thresholds
  ThresholdRatio = .4;          % Low thresh is this fraction of the high.
  
  % Design the filters - a gaussian and its derivative
  
  pw = 1:30; % possible widths
  ssq = sigma*sigma;
  width = max(find(exp(-(pw.*pw)/(2*sigma*sigma))>GaussianDieOff));
  if isempty(width)
    width = 1;  % the user entered a really small sigma
  end

  t = (-width:width);
  gau = exp(-(t.*t)/(2*ssq))/(2*pi*ssq);     % the gaussian 1D filter

  % Find the directional derivative of 2D Gaussian (along X-axis)
  % Since the result is symmetric along X, we can get the derivative along
  % Y-axis simply by transposing the result for X direction.
  [x,y]=meshgrid(-width:width,-width:width);
  dgau2D=-x.*exp(-(x.*x+y.*y)/(2*ssq))/(pi*ssq);
  
  % Convolve the filters with the image in each direction
  % The canny edge detector first requires convolution with
  % 2D gaussian, and then with the derivitave of a gaussian.
  % Since gaussian filter is separable, for smoothing, we can use 
  % two 1D convolutions in order to achieve the effect of convolving
  % with 2D Gaussian.  We convolve along rows and then columns.

  %smooth the image out
  aSmooth=imfilter(a,gau,'conv','replicate');   % run the filter accross rows
  aSmooth=imfilter(aSmooth,gau','conv','replicate'); % and then accross columns
  
  %apply directional derivatives
  ax = imfilter(aSmooth, dgau2D, 'conv','replicate');
  ay = imfilter(aSmooth, dgau2D', 'conv','replicate');

  mag = sqrt((ax.*ax) + (ay.*ay));
  magmax = max(mag(:));
  if magmax>0
    mag = mag / magmax;   % normalize
  end
  
%   figure;
%   %rotular2((1-mag)*4095);
%   imshow(gau,[]);
  
  return
  
  [counts,x]=imhist(mag, 64);
    highThresh = min(find(cumsum(counts) > PercentOfPixelsNotEdges*m*n)) / 64;
    
    %[highThresh imagen_filtrada]=bivaluada(mag*4095,1);
    highThresh = highThresh/255;
    lowThresh = ThresholdRatio*highThresh;
    
    thresh = [lowThresh highThresh];
    
    
  idxStrong = [];  
  for dir = 1:4
    idxLocalMax = cannyFindLocalMaxima(dir,ax,ay,mag);
    idxWeak = idxLocalMax(mag(idxLocalMax) > lowThresh);
    e(idxWeak)=1;
    idxStrong = [idxStrong; idxWeak(mag(idxWeak) > highThresh)];
  end
  
  rstrong = rem(idxStrong-1, m)+1;
  cstrong = floor((idxStrong-1)/m)+1;
  e = bwselect(e, cstrong, rstrong, 8);
  e = bwmorph(e, 'thin', 1);
  
  
  
  
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Local Function : cannyFindLocalMaxima
%
function idxLocalMax = cannyFindLocalMaxima(direction,ix,iy,mag);
%
% This sub-function helps with the non-maximum supression in the Canny
% edge detector.  The input parameters are:
% 
%   direction - the index of which direction the gradient is pointing, 
%               read from the diagram below. direction is 1, 2, 3, or 4.
%   ix        - input image filtered by derivative of gaussian along x 
%   iy        - input image filtered by derivative of gaussian along y
%   mag       - the gradient magnitude image
%
%    there are 4 cases:
%
%                         The X marks the pixel in question, and each
%         3     2         of the quadrants for the gradient vector
%       O----0----0       fall into two cases, divided by the 45 
%     4 |         | 1     degree line.  In one case the gradient
%       |         |       vector is more horizontal, and in the other
%       O    X    O       it is more vertical.  There are eight 
%       |         |       divisions, but for the non-maximum supression  
%    (1)|         |(4)    we are only worried about 4 of them since we 
%       O----O----O       use symmetric points about the center pixel.
%        (2)   (3)        


[m,n,o] = size(mag);

% Find the indices of all points whose gradient (specified by the 
% vector (ix,iy)) is going in the direction we're looking at.  

switch direction
 case 1
  idx = find((iy<=0 & ix>-iy)  | (iy>=0 & ix<-iy));
 case 2
  idx = find((ix>0 & -iy>=ix)  | (ix<0 & -iy<=ix));
 case 3
  idx = find((ix<=0 & ix>iy) | (ix>=0 & ix<iy));
 case 4
  idx = find((iy<0 & ix<=iy) | (iy>0 & ix>=iy));
end

% Exclude the exterior pixels
if ~isempty(idx)
  v = mod(idx,m);
  extIdx = find(v==1 | v==0 | idx<=m | (idx>(n-1)*m));
  idx(extIdx) = [];
end

ixv = ix(idx);  
iyv = iy(idx);   
gradmag = mag(idx);

% Do the linear interpolations for the interior pixels
switch direction
 case 1
  d = abs(iyv./ixv);
  gradmag1 = mag(idx+m).*(1-d) + mag(idx+m-1).*d; 
  gradmag2 = mag(idx-m).*(1-d) + mag(idx-m+1).*d; 
 case 2
  d = abs(ixv./iyv);
  gradmag1 = mag(idx-1).*(1-d) + mag(idx+m-1).*d; 
  gradmag2 = mag(idx+1).*(1-d) + mag(idx-m+1).*d; 
 case 3
  d = abs(ixv./iyv);
  gradmag1 = mag(idx-1).*(1-d) + mag(idx-m-1).*d; 
  gradmag2 = mag(idx+1).*(1-d) + mag(idx+m+1).*d; 
 case 4
  d = abs(iyv./ixv);
  gradmag1 = mag(idx-m).*(1-d) + mag(idx-m-1).*d; 
  gradmag2 = mag(idx+m).*(1-d) + mag(idx+m+1).*d; 
end
idxLocalMax = idx(gradmag>=gradmag1 & gradmag>=gradmag2); 