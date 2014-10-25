function uh = hermiteslow(n, X, U, dU, x)
%HERMITESLOW - Calcule et évalue une interpolation cubique par morceaux.
%
%   Cette fonction calcule les coefficients des n polynomes d'interpolation
%   cubique P(i) selon la règle suivante :
%   P(i)  = U(i) en X(i) et U(i+1) en X(i+1)
%   P'(i) = dU(i) en X(i) et dU(i+1) en X(i+1),
%   avec X, U et dU de taille n+1.
%   Elle évalue ensuite chaque polynome sur une section du vecteur x,
%   correspondant à l'intervalle sur lequel le polynome est défini.
%
%   HERMITE(n, X, U, dU, x)
%       n  = le nombre d'intervalles
%       X  = les n+1 abscisses d'interpolation, délimitant les n intervalles
%       U  = les valeurs de la fonction à interpoler en chaque X
%       dU = les valeurs de la dérivée première de la fonction à interpoler
%               en chaque X
%       x  = le vecteur d'abscisses sur lequel évaluer les polynomes
%               d'interpolation

uh = zeros(length(x));

for i=1:n
    % On résoud un système pour trouver les 4 coefficients du i-ème
    % polynome.
    % Ce système correspond aux règles posées :
    % P(i)  = U(i) en X(i) et U(i+1) en X(i+1)
	% P'(i) = dU(i) en X(i) et dU(i+1) en X(i+1).
    A = [   X(i)^3,     X(i)^2,   X(i),   1 ;
          3*X(i)^2,   2*X(i),     1,      0 ;
            X(i+1)^3,   X(i+1)^2, X(i+1), 1 ;
          3*X(i+1)^2, 2*X(i+1),   1,      0 ];
    b = [ U(i) ; dU(i) ; U(i+1) ; dU(i+1) ];
    
    p = A\b;
    
    % Ensuite, on évalue ce polynome sur les abscisses du vecteur x qui
    % sont comprises entre X(i) et X(i+1).
    for j=1:length(x)
        if X(i) <= x(j) && x(j) <= X(i+1)
            uh(j) = poly3eval(p, x(j));
        end
    end
end

% Une boucle sur n qui imbrique une boucle sur m=length(x) nous donne une
% complexité O(n*m).
% On pourrait réduire celle-ci en supposant que les abscisses de x sont
% nécessairement ordonnées, ou en les ordonnant au départ.
% Ceci est fait dans "hermite.m".

end

function uh = poly3eval(p, x)
    uh = p(1)*x^3 + p(2)*x^2 + p(3)*x + p(4);
end