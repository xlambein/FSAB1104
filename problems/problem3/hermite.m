function uh = hermite(n, X, U, dU, x)
%HERMITE - Calcule et evalue une interpolation cubique par morceaux.
%
%   Cette fonction calcule les coefficients des n polynomes d'interpolation
%   cubique P(i) selon la regle suivante :
%   P(i)  = U(i) en X(i) et U(i+1) en X(i+1)
%   P'(i) = dU(i) en X(i) et dU(i+1) en X(i+1),
%   avec X, U et dU de taille n+1.
%   Elle evalue ensuite chaque polynome sur une section du vecteur x,
%   correspondant a l'intervalle sur lequel le polynome est defini.
%
%   HERMITE(n, X, U, dU, x)
%       n  = le nombre d'intervalles
%       X  = les n+1 abscisses d'interpolation, delimitant les n intervalles
%       U  = les valeurs de la fonction a interpoler en chaque X
%       dU = les valeurs de la derivee premiere de la fonction a interpoler
%               en chaque X
%       x  = le vecteur d'abscisses sur lequel evaluer les polynomes
%               d'interpolation

uh = zeros(length(x));
p = zeros(4, n);

for i=1:n
    % On resoud un systeme pour trouver les 4 coefficients du i-eme
    % polynome.
    % Ce systeme correspond aux regles posees :
    % P(i)  = U(i) en X(i) et U(i+1) en X(i+1)
	% P'(i) = dU(i) en X(i) et dU(i+1) en X(i+1).
    A = [   X(i)^3,     X(i)^2,   X(i),   1 ;
          3*X(i)^2,   2*X(i),     1,      0 ;
            X(i+1)^3,   X(i+1)^2, X(i+1), 1 ;
          3*X(i+1)^2, 2*X(i+1),   1,      0 ];
    b = [ U(i) ; dU(i) ; U(i+1) ; dU(i+1) ];
    
    p(:, i) = A\b;
end

% Dans une premiere version de hermite, l'evaluation etait faite sur le
% vecteur x a l'interieur d'une boucle sur n :
% for j=1:length(x)
%     if X(i) <= x(j) && x(j) <= X(i+1)
%         uh(j) = poly3eval(p, x(j));
%     end
% end
% Une telle implementation nous donne une complexite de O(n*m), avec
% m=length(x).
% On peut reduire celle-ci en supposant que les abscisses de x sont
% necessairement ordonnees, ou en les ordonnant au depart. C'est ce que
% j'ai decide de faire ici :

x = sort(x);

% Ensuite, il suffit d'evaluer en chaque x(i) le j-eme polynome, en
% incrementant j des que x(i) depasse X(j+1) :
j = 1;
for i=1:length(x)
    if x(i) > X(j+1)
        j = j+1;
    end
    uh(i) = poly3eval(p(:, j), x(i));
end

% On obtient une complexite de n pour le calcul des coefficients, m*log(m)
% pour sort(x) et m pour l'evaluation des polynomes sur le vecteur x.
% La complexite de hermite est donc O(n + m*log(m)).

end

function uh = poly3eval(p, x)
    uh = p(1)*x^3 + p(2)*x^2 + p(3)*x + p(4);
end
