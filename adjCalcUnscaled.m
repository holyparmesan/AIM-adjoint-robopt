function adj_maps = adjCalc(E12, E21, index, varargin)
% The first three inputs are structs, to be taken directly from Lumerical's varFDTD.
% The remaining inputs are contours.

E12dat = unFlattenData(E12.E,E12.x,E12.y);
E21dat = unFlattenData(E21.E,E21.x,E21.y);
idx_x = unFlattenData(index.index_x,index.x,index.y);
idx_y = unFlattenData(index.index_y,index.x,index.y);

D12dat = zeros(size(E12dat));
D12dat(:,:,1,:) = E12dat(:,:,1,:) .* (idx_x.^2);
D12dat(:,:,2,:) = E12dat(:,:,2,:) .* (idx_y.^2);
D21dat = zeros(size(E21dat));
D21dat(:,:,1,:) = E21dat(:,:,1,:) .* (idx_x.^2);
D21dat(:,:,2,:) = E21dat(:,:,2,:) .* (idx_y.^2);

n_hi = max(idx_x(:));
n_lo = min(idx_x(:));
deps = n_hi^2 - n_lo^2;
dieps = n_hi^(-2) - n_lo^(-2);

[xsm, ysm] = meshgrid(E12.x, E12.y);

adj_maps = cell(nargin-3,1);

for i = 1:nargin-3
    contour = varargin{i};
    N = size(contour, 1);
    Nlambs = size(E12dat,4);
    adj_maps{i} = zeros(N,Nlambs);
    
    E12paras = zeros(N,Nlambs);
    E21paras = zeros(N,Nlambs);
    D12perps = zeros(N,Nlambs);
    D21perps = zeros(N,Nlambs);

    deltas = zeros(N,2);
    deltas(2:N-1,:) = contour(3:N,:) - contour(1:N-2,:);
    deltas(1,:) = contour(2,:) - contour(1,:);
    deltas(N,:) = contour(N,:) - contour(N-1,:);
    thetas = atan2(deltas(:,2),deltas(:,1));
    
    % approach: interpolate all four fields, then multiply
    for j = 1:N
        
        rxi = find(E12.x > contour(j,1),1); % rough x index
        ryi = find(E12.y > contour(j,2),1); % rough y index
                
        % yeah i know, i'm sorry :'(
        for k = 1:Nlambs
            E12paras(j,k) = interp2(xsm(ryi-4:ryi+3,rxi-4:rxi+3), ysm(ryi-4:ryi+3,rxi-4:rxi+3), (E12dat(rxi-4:rxi+3,ryi-4:ryi+3,1,k)*cos(thetas(j)) + E12dat(rxi-4:rxi+3,ryi-4:ryi+3,2,k)*sin(thetas(j))).', contour(j,1), contour(j,2));
            E21paras(j,k) = interp2(xsm(ryi-4:ryi+3,rxi-4:rxi+3), ysm(ryi-4:ryi+3,rxi-4:rxi+3), (E21dat(rxi-4:rxi+3,ryi-4:ryi+3,1,k)*cos(thetas(j)) + E21dat(rxi-4:rxi+3,ryi-4:ryi+3,2,k)*sin(thetas(j))).', contour(j,1), contour(j,2));
            D12perps(j,k) = interp2(xsm(ryi-4:ryi+3,rxi-4:rxi+3), ysm(ryi-4:ryi+3,rxi-4:rxi+3), (D12dat(rxi-4:rxi+3,ryi-4:ryi+3,1,k)*sin(thetas(j)) - D12dat(rxi-4:rxi+3,ryi-4:ryi+3,2,k)*cos(thetas(j))).', contour(j,1), contour(j,2));
            D21perps(j,k) = interp2(xsm(ryi-4:ryi+3,rxi-4:rxi+3), ysm(ryi-4:ryi+3,rxi-4:rxi+3), (D21dat(rxi-4:rxi+3,ryi-4:ryi+3,1,k)*sin(thetas(j)) - D21dat(rxi-4:rxi+3,ryi-4:ryi+3,2,k)*cos(thetas(j))).', contour(j,1), contour(j,2));
        end
    end
    
    adj_maps{i} = E12paras.*E21paras*deps - D12perps.*D21perps*dieps;
    % This still needs to be multiplied by 2pi*i*f*e0 / N!
    
end

end