function [uh] = lagrange(x,X,U)

n = min(length(X),length(U));  % Evite un message d'erreur :-)
phi = ones(n,length(x));       % Preallocation (important !)
for i =1:n
    for j = [1:(i-1) (i+1):n]
          phi(i,:) = phi(i,:).*(x-X(j))/(X(i)-X(j));
    end
end
uh = U * phi;                  % Produit matriciel (plus rapide !)


end

