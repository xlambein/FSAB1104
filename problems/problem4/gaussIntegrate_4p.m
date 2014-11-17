function [ Ih ] = gaussIntegrate( L, z, error )
%GAUSSINTEGRATE Summary of this function goes here
%   Detailed explanation goes here

x = [-sqrt(3/7 + 2/7*sqrt(6/5)),...
     -sqrt(3/7 - 2/7*sqrt(6/5)),...
      sqrt(3/7 - 2/7*sqrt(6/5)),...
      sqrt(3/7 + 2/7*sqrt(6/5))];

w = [(18 - sqrt(30))/36,...
     (18 + sqrt(30))/36,...
     (18 + sqrt(30))/36,...
     (18 - sqrt(30))/36];

W = w'*w;
X = ones(4, 1)*x;
Y = x'*ones(1, 4);

U = z(X, Y);

Ih = L*L*sum(sum(W*U));

end

