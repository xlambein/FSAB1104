% Exercice 9 (APE2)

function exo9
syms x

estimateError('cosh', (-1 + (0:4)/2))
estimateError('sinh', (-1 + (0:4)/2))
estimateError(matlabFunction(cos(x) + sin(x)), (-pi/2 + pi*(0:4)/4))

end

function err = estimateError(f, X)
% ESTIMATEERROR  Estimates error between a function and its polynomial
%                approximation using points (X, f(X)).

x = linspace(X(1), X(end), 1000);
U = feval(f, X);
u = feval(f, x);

uh = polyval(polyfit(X, U, length(X)-1), x);

err = max(u-uh);
end