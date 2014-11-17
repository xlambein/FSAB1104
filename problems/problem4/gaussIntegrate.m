function [ Ih ] = gaussIntegrate( L, z, error )
%GAUSSINTEGRATE Integration with Gauss-Legendre over a surface.
%   GAUSSINTEGRATE(L, z, e) is the integral of the two-variables function
%   z over the surface [-L, L]x[-L, L] with an error smaller than e.

% We ensure that the values make sense
L = abs(L(1));
error = abs(error(1));

% To compute the best value for Ih, we proceed recursively : at each
% iteration, we compute a more precise integral with Gauss-Legendre, by
% dividing the surface of integration by two by four (here done by
% multiplying n by two). Then, we obtain an even better value by using
% a Richardson extrapolation with the previously computed values.
% If the difference between that value and the one from the previous
% iteration is bigger than the error, we repeat. Otherwise, it means enough
% iterations have been done and thus we just return that best value we just
% computed.

% First, we need to calculate the lowest-precision approximation. This is
% done by using the Gauss-Legendre method over the full domain
% [-L, L]x[-L, L].
I00 = gl(L, 1, z);

% Then, we call the recursive function :
Ih = precisionStep(L, 2, z, error, I00);

end

function Ih = precisionStep(L, n, z, error, prevI)

% Here we compute the integral with Gauss-Legendre. 
Ii0 = gl(L, n, z);

% Then we compute the next line of the Richardson extrapolation. The last
% element of this vector is our current best value for Ih.
Ii = nextRichardson(prevI, Ii0);

% If the difference between the best result and the previous best result is
% still bigger than error, then it means we need to do (at least) another
% iteration.
if abs(Ii(end)-prevI(end)) > error
    Ih = precisionStep(L, 2*n, z, error, Ii);
else
    Ih = Ii(end);
end

end

function I = nextRichardson(prevI, I0)
% This function is used to compute the next vector of values used in
% the Richardson extrapolation, using the previous one, as well as the
% first value of the next vector.

I = zeros(length(prevI)+1, length(I0));
I(1) = I0;

for j=2:size(I)
    a = 2^(4+j);
    I(j) = (a*I(j-1) - prevI(j-1))/(a - 1);
end

end

function Ih = gl(L, n, z)

h = L/n;
s3 = sqrt(3);

% The integration is done by evaluating the function z at four points in
% n^2 surfaces of side h=L/n. These four points are, of course, the points
% given by the Gauss-Legendre formula : (+-h/sqrt(3), +-h/sqrt(3)).

c = linspace(-L+h, L-h, n);
X = repmat([c - h/s3, c + h/s3], 2*n, 1);
Ih = sum(sum(z(X, X')))*h*h;

end