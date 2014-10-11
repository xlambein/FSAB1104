function exo26

end

function B = Bspline(T, t, p, i)

if p==0
    B = 

B = (t - T(i))/(T(i+p) - T(i))*Bspline(T, t, p-1, i)...
  + (T(i+p+1) - t)/(T(i+p+1) - T(i+1))*Bspline(T, t, p-1, i+1);

end