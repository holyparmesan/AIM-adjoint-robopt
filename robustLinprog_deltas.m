function bestdelta = robustLinprog_deltas(f0s, grads, inflim, deltlim)
% f0s is k x 1, grads is k x N, inflim is 1 x 1.
%
% We have an N-dimensional variable, x
% We have k linear functions, the ith of which is f_i(x) = f0s(i) + grads(i,:) * x
% We seek to find the x which minimizes max_i f_i(x)
% subject to |x|_inf <= inflim
%
% This can be re-cast as an N+1-dimensional optimization problem in [t x]
% where we seek to minimize t
% subject to t >= f_i(x) for all i
% and also -deltlim <= x_j - x_j+1 <= deltlim for all j
% and also -inflim <= x_j <= inflim for all j

[k, N] = size(grads);

f = [1; zeros(N,1)];
A = [-1*ones(k,1), grads];
b = -1*f0s;
if numel(inflim) == 1
    lim = [Inf;inflim*ones(N,1)];
else
    lim = [Inf;inflim];
end

if numel(deltlim) == 1
    dlim = deltlim*ones(2*(N-1),1);
else
    assert(numel(deltlim) == N-1);
    dlim = [deltlim;deltlim];
end

dmat = [zeros(N-1,1) eye(N-1) zeros(N-1,1)] - [zeros(N-1,2) eye(N-1)];

% minimize f.' * x s.t. A*x <= b

soln = linprog(f,[A;dmat;-dmat],[b;dlim],[],[],-lim,lim);
bestdelta = soln(2:end);

end

