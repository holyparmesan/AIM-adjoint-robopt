function [delts, norms] = delts_norms(varargin)

delts = cell(nargin,1);
norms = cell(nargin,1);
for k = 1:nargin
    [delts{k}, norms{k}] = delt_norm(varargin{k});
end