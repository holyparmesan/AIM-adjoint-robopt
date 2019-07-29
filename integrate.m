function sum = integrate(vector, contour)

sum = 0;

N = size(contour,1);
deltas = zeros(N,2);
deltas(2:N-1,:) = contour(3:N,:) - contour(1:N-2,:);
deltas(1,:) = contour(2,:) - contour(1,:);
deltas(N,:) = contour(N,:) - contour(N-1,:);

lengths = vecnorm(deltas.').' / 2;

for i = 1:N
    sum = sum + lengths(i)*vector(i); % trapezoiding as HECK
end

end

