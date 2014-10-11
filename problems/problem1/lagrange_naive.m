function [uh] = lagrange_naive(x,X,U)

n = min(length(X),length(U));
for j = 1:length(x)
   uh(j) = 0;
   for i =1:n 
      phi(i) = 1;
      for k =1:n
         if (i ~= k)
            phi(i) = phi(i)*(x(j)-X(k))/(X(i)-X(k));
         end
      end
      uh(j) = uh(j) + U(i) * phi(i);
   end
end

end
