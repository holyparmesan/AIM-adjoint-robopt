function adj_maps = adjCalc(xs, ys, E12dat, E21dat, idx_x, idx_y, f, N_mode, distanceScale, varargin)
% The next two inputs are E-field maps with dimensions {xs, ys, 1, freqs, 3}
% The next two are effective index maps with dimensions {xs, ys, 1, freqs}
% f and N_mode are for scaling purposes, and should be 1D
% distanceScale is a boolean
% the rest are contours

%E12dat = unFlattenData(E12.E,E12.x,E12.y);
%E21dat = unFlattenData(E21.E,E21.x,E21.y);
%idx_x = unFlattenData(index.index_x,index.x,index.y);
%idx_y = unFlattenData(index.index_y,index.x,index.y);

D12dat = zeros(size(E12dat));
D12dat(:,:,:,:,1) = E12dat(:,:,:,:,1) .* (idx_x.^2);
D12dat(:,:,:,:,2) = E12dat(:,:,:,:,2) .* (idx_y.^2);
D21dat = zeros(size(E21dat));
D21dat(:,:,:,:,1) = E21dat(:,:,:,:,1) .* (idx_x.^2);
D21dat(:,:,:,:,2) = E21dat(:,:,:,:,2) .* (idx_y.^2);

n_hi = max(idx_x(:));
n_lo = min(idx_x(:));
deps = n_hi^2 - n_lo^2;
dieps = n_hi^(-2) - n_lo^(-2);

[xsm, ysm] = meshgrid(xs, ys);

adj_maps = cell(nargin-9,1);

for i = 1:nargin-9
    contour = varargin{i};
    N = size(contour, 1);
    Nlambs = size(E12dat,4);
    adj_maps{i} = zeros(N,Nlambs);
    
    E12paras = zeros(N,Nlambs);
    E21paras = zeros(N,Nlambs);
    D12perps = zeros(N,Nlambs);
    D21perps = zeros(N,Nlambs);

    [deltas, ~] = delt_norm(contour);
    thetas = atan2(deltas(:,2),deltas(:,1));
    
    % approach: interpolate all four fields, then multiply
    for j = 1:N
        
        rxi = find(xs > contour(j,1),1); % rough x index
        ryi = find(ys > contour(j,2),1); % rough y index
                
        % yeah i know, i'm sorry :'(
        for k = 1:Nlambs
            E12paras(j,k) = interp2(xsm(ryi-3:ryi+2,rxi-3:rxi+2), ysm(ryi-3:ryi+2,rxi-3:rxi+2), (E12dat(rxi-3:rxi+2,ryi-3:ryi+2,1,k,1)*cos(thetas(j)) + E12dat(rxi-3:rxi+2,ryi-3:ryi+2,1,k,2)*sin(thetas(j))).', contour(j,1), contour(j,2));
            E21paras(j,k) = interp2(xsm(ryi-3:ryi+2,rxi-3:rxi+2), ysm(ryi-3:ryi+2,rxi-3:rxi+2), (E21dat(rxi-3:rxi+2,ryi-3:ryi+2,1,k,1)*cos(thetas(j)) + E21dat(rxi-3:rxi+2,ryi-3:ryi+2,1,k,2)*sin(thetas(j))).', contour(j,1), contour(j,2));
            D12perps(j,k) = interp2(xsm(ryi-3:ryi+2,rxi-3:rxi+2), ysm(ryi-3:ryi+2,rxi-3:rxi+2), (D12dat(rxi-3:rxi+2,ryi-3:ryi+2,1,k,1)*sin(thetas(j)) - D12dat(rxi-3:rxi+2,ryi-3:ryi+2,1,k,2)*cos(thetas(j))).', contour(j,1), contour(j,2));
            D21perps(j,k) = interp2(xsm(ryi-3:ryi+2,rxi-3:rxi+2), ysm(ryi-3:ryi+2,rxi-3:rxi+2), (D21dat(rxi-3:rxi+2,ryi-3:ryi+2,1,k,1)*sin(thetas(j)) - D21dat(rxi-3:rxi+2,ryi-3:ryi+2,1,k,2)*cos(thetas(j))).', contour(j,1), contour(j,2));
        end
    end
    
    % we need to multiply by 2pi*i*f*e0 / 4N
    % and we need to integrate along the contour length
    % each point contributes an "effective length" of deltas() / 2
    if (distanceScale)
        scaling = vecnorm(deltas.').' * (f./N_mode).' * (pi * 8.8142e-12j/ 4); 
    else
        scaling = ones(N,1) * (f./N_mode).' * (2*pi * 8.8142e-12j/ 4);
    end
    adj_maps{i} = (E12paras.*E21paras*deps - D12perps.*D21perps*dieps) .* scaling;
    
end

end

function res = fitinterp2(xs, ys, data, x, y)
    params = [ones(numel(xs),1) xs(:) ys(:)];
    slopes = params \ data(:);
    res = [1 x y] * slopes;
end