% Exercice 9 (APE2)

function exo9
syms x
f3 = cos(x) + sin(x);

estimateError('cosh', (-1 + (0:4)/2))
estimateError('sinh', (-1 + (0:4)/2))
estimateError(f3, (-pi/2 + pi*(0:4)/4))

end

function err = estimateError(f, X)
U = feval(f, X);

x = linspace(X(1), X(end), 1000);
u = feval(f, x);
uh = polyval(polyfit(X, U, length(X)-1), x);

err = max(u-uh);
end