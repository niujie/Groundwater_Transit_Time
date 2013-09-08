function [d, t, p] = p2ldist(x, y)
% P2LDIST Distance from points to polygonal line
% [D, T, P] = P2LDIST (X, A, B) returns the distance from the point X to 
% the polygonal line defined by array Y. Rows correspond to points. X must be a 
% row vector. X and Y must have the same number of columns. Return arguments 
% are column vectors D and T, which hold the squared distances and the normalized
% parameter values achieve them. Array P contains the closest point on Y to X.

a = y(1:(end - 1), :);
b = y(2:end, :);
[d, t, p] = p2sdist(x, a, b);
return

function [d, t, p] = p2sdist(x, a, b)
% P2SDIST Distance from points to segments
% [D, T, P] = P2SDIST(X, A, B) returns the distance from the points in X to 
% the segments with endpoints A and B. Both A and B must be the same size.
% Rows of correspond to points. X can have either a single row or as many rows 
% as A and B. Return arguments are column vectors D and T that contain, 
% respectively, the squared distances and the normalized parameter values that 
% achieve them. Array P contains the closest point on each segment AB to each 
% point X.

v = b - a;
z = bsxfun(@minus, x, a);
u = sum(v.^2, 2);
t = sum(bsxfun(@times, z, v), 2);
i = u > eps;
t(i) = t(i)./u(i);
t(~i) = 0.0;
t = max(0.0, min(t, 1.0));
p = bsxfun(@times, t, v);
z = z - p;
d = sum(z.^2, 2);
p = bsxfun(@plus, a, p);
return