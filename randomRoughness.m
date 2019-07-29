function [flucts, corrmat] = randomRoughness(contour, CL, k)
% contour is an Nx2 array of points
% flucts is an Nxk array of fluctuations
% flucts is multivariate Gaussian; correlation goes as exp(-(dist / CL)^2)

xdiffs = contour(:,1) - contour(:,1).';
ydiffs = contour(:,2) - contour(:,2).';

corrmat = exp(-(xdiffs.^2 + ydiffs.^2)/CL^2);

flucts = mvnrnd(zeros(1,size(corrmat,1)),corrmat,k);

end

