function [ Ih ] = gaussIntegrate( L, z, error )
%GAUSSINTEGRATE Summary of this function goes here
%   Detailed explanation goes here

L = abs(L);
error = abs(error);

% To calculate the best value for Ih, we proceed by using Romberg's method
% recursively, as many time as needed for each sub-intervals, until the
% error estimated in each sub-interval is less than ERROR.

% First, we need to calculate the lowest-precision approximation. This is
% done by using the Gauss-Legendre method over the domain [-L, L]x[-L, L],
% centered in (0, 0).
I00 = gl(L, z, 0, 0);

% Then, we call the recursive function with a step divided by two :
Ih = precisionStep(L/2, z, error/4, I00, 0, 0);

end

function I = precisionStep(h, z, error, prevI, Dx, Dy)

I = 0;

% We consider a square of side 2*h, centered in (Dx, Dy), whose current
% best integral is contained in prevI(end).
% We are going to subdivise this square in 4 squares of side h, centered
% relatively to this square in (-h, -h), (-h, h), (h, -h), (h, h). Then, we
% are going to integrate these 4 squares using Gauss-Legendre. These
% integrals, along with the integrals of the bigger square contained in
% prevI, will allow us to calculate more accurate integrals with Romberg's
% method.

% These are the coordinates of the centers of the 4 little squares :
x = [Dx-h, Dx-h, Dx+h, Dx+h];
y = [Dy-h, Dy+h, Dy-h, Dy+h];

% Here we integrate the 4 squares using Gauss-Legendre.
for k=1:4
    Ii0(k) = gl(h, z, x(k), y(k));
end

% This function calculates a better approximation using Romberg's method
% with the integrals we just calculated, and the different precision
% integrals of the bigger square stored in prevI.
Ii = nextRichardson(prevI, Ii0);

% Finally, we check the size of the error for each little square.
for k=1:4
    if abs(Ii(end, k)-Ii(end-1, k)) > error
        % If the difference between the best approximation and the second
        % best approximation is bigger than the error, it means we need to
        % do Romberg's method with a lower step h/2.
        % Note that we divide the error by 4 : this is because, since the
        % final integral is the sum of the 4 smaller integral, the error in
        % the final integral is equal to 4 times the biggest error in the
        % smaller integrals. Thus, we need the errors in these smaller
        % integrals to be 4 times smaller.
        % For the same kind of reason, we divide Ii(:, k) by 4.
        % Once calculated, that value is added to the total integral of the
        % big square
        Istar = precisionStep(h/2, z, error/4, Ii(:, k)./4, x(k), y(k));
        I = I + Istar;
    else
        %length(Ii(:, k))
        % If the error is reasonable, then we simply add the best value of
        % the integral of the little square to the total integral of the
        % big square.
        I = I + Ii(end, k);
    end
end

end

function I = nextRichardson(prevI, I0)

I = zeros(length(prevI)+1, length(I0));
I(1, :) = I0;

for j=2:size(I)
    a = 4^(4+j);
    I(j, :) = (a*I(j-1, :) - prevI(j-1))/(a - 1);
end

end

function Ih = gl(h, z, Dx, Dy)

s3 = 1/sqrt(3);
X = [-s3 -s3 s3 s3];
Y = [-s3 s3 -s3 s3];
W = [1 1 1 1];

U = z(X*h+Dx, Y*h+Dy);
Ih = sum(W.*U)*h*h;

end