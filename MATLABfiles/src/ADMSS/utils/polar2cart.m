function cartIm = polar2cart(rtheta,cartSize,bkgColor,deltaCat)
% Forms the cartesian form of the IVUS image
% The polar image is without catether

[x,y] = size(rtheta);

% minimun and maximun depth in pixels
min_depth = deltaCat;
max_depth = cartSize/2; 

% delta of the radius
dr = (max_depth - min_depth)/x;
r = [min_depth:dr:max_depth-dr];

% create a linear space from 0 to 2pi angles and y+1 lines
t = linspace(0,2*pi,y+1);
t = t(1:end-1);

% create a mesh in order to convert it into the polar form
[theta,rho] = meshgrid(t,r);
[x,y] = pol2cart(theta,rho);

% convert into polar form
ti = -(cartSize/2):((cartSize/2)-1);

% make the interpolation
[XI,YI] = meshgrid(ti,ti);
ZI = griddata(y,x,rtheta,XI,YI);
%Codigo mio>
Mask=(XI.^2 + YI.^2) <= (deltaCat.^2);
ZI(Mask==1)=Inf;
% filling the background
tf = isnan(ZI);
tf2 = isinf(ZI);
h = find(tf==1);
h2 = find(tf2==1);
cartIm = ZI;
cartIm(h) = bkgColor;
cartIm(h2) = bkgColor;
cartIm = fliplr(cartIm);
cartIm = flipud(cartIm); % to delete to match with dicom

