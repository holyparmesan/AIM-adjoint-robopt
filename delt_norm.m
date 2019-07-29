function [delts, norms] = delt_norm(cont)

delts = zeros(size(cont));
delts(1,:) = cont(2,:) - cont(1,:);
delts(2:end-1,:) = cont(3:end,:) - cont(1:end-2,:);
delts(end,:) = cont(end,:) - cont(end-1,:);
   
norms_raw = [delts(:,2), -1*delts(:,1)]; 
norms = norms_raw ./ vecnorm(norms_raw.').';