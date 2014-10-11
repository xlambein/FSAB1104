% Exercice 12 (APE2)
function t = exo12(T, x, y)

X = [-1 -1 1 1];
Y = [-1 1 1 -1];

C = times(X, Y);

M = [ones(1, 4) ; X(1:end) ; Y(1:end) ; C(1:end)]';

a = M\T';

t = a(1) + x*a(2) + y*a(3) + x.*y*a(4);

% BONUS: graphe 3D de l'interpolation
% [x, y] = meshgrid(-1:.1:1, -1:.1:1);
% th = a(1) + x*a(2) + y*a(3) + x.*y*a(4);
% surf(x, y, th)

end