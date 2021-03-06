function [uh] = lagrange_super(x,X,U)
%LAGRANGE_SUPER - Interpolation polynomiale par polynomes de Lagrange et
%évaluation sur un vecteur
%
%   Cette fonction calcule les coefficients d'un polynome de degré
%   min(length(X),length(U)) qui interpole les points (X(i), U(i)), puis
%   évalue ce polynome sur le vecteur x.
%
%   uh = LAGRANGE_SUPER(x,X,U)
%
%   See also LAGRANGE, LAGRANGE_NAIVE, POLYFIT, POLYVAL.

n = min(length(X),length(U));
m = length(x);

% Calcul de la partie rouge
rouge = ones(1, m);
for i = 1:n
    rouge = rouge.*(x - X(i)); % O(m*n)
end

% Calcul de uh
uh = zeros(1, m);
for i = 1:n
    % Calcul du i-ème terme en bleu
    bleu = U(i);
    for j = [1:(i-1) (i+1):n]
        bleu = bleu/(X(i) - X(j)); % O(n^2)
    end
    
    % C'est ici qu'un problème pourrait survenir, si une des valeurs de x
    % était exactement égale à X(i). En effet, matlab travaillant avec des
    % nombres et non des symboles, il n'est pas possible d'éviter une
    % division par zéro en annulant le numérateur et le dénominateur.
    %
    % La façon dont on adresse ce problème est la suivante :
    % au lieu d'utiliser linspace normalement, on se débrouille pour que
    % les bornes de x ne soient pas exactement égales aux bornes de X, mais
    % à la place diffèrent d'une valeur très petite, qui, elle, peut être
    % annulée par son inverse. Cette valeur nous est fournit par matlab à
    % l'aide de la fonction 'eps'. Ainsi, en utilisant linspace de la façon
    % suivante :
    % 'linspace(eps, 5*(1-eps), m)',
    % on s'assure qu'il n'y ait pas de division par zéro.
    % Cependant, tous les problèmes ne sont pas évités : en procédant
    % ainsi, le programme est peu robuste et rapidement victime de pertes
    % de précision dues aux nombres flottants. C'est pour cette raison que,
    % pour des valeurs de n > 60, les bords de l'interpolation deviennent
    % progressivement de plus en plus éloignés de la fonction interpolée.
    
    uh = uh + bleu./(x - X(i)); % O(m*n)
end
uh = uh .* rouge; % O(m)

% La complexité des opérations importantes est indiquée en commentaire à
% côté de celles-ci. Au total, on en note trois différentes :
% O(m), O(n^2) et O(m*n).
% Avec m > n, la borne supérieur de complexité de l'algorithme est O(m*n),
% ce qui est bien meilleur que O(m*n^2) dans lagrange.m

end

