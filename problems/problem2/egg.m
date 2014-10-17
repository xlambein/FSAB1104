function [x, y, z] = egg(top,bottom,dt,mode)

s2 = sqrt(2)/2;

T = [0 0 0 1 1 2 2 3]; %nt = 7
S = [0 0 0 1 1 2 2 3 3 4 4 5]; %ns = 11
R = [0 1 1 1 0];
H = [-1*bottom -1*bottom 0 1*top 1*top];


p = 2;
Xc = [1 1 0 -1 -1 -1 0 1 1] ;
Yc = [0 1 1 1 0 -1 -1 -1 0] ;
Zc = ones(size(Xc));
Wc = [1 s2 1 s2 1 s2 1 s2 1];
X = Xc' * R;
Y = Yc' * R;
Z = Zc' * H;
W = Wc' * [1 s2 1 s2 1];

nt = length(T) - 1;
t = [T(p+1):dt:T(nt-p+1)];
for i=0:nt-p-1
  Bt(i+1,:) = b(t,T,i,p);
end

ns = length(S) - 1;
s = [S(p+1):dt:S(ns-p+1)];
for i=0:ns-p-1
  Bs(i+1,:) = b(s,S,i,p);
end

w = Bs' * W * Bt;
x = Bs' * (W .* X) * Bt ./ w;
y = Bs' * (W .* Y) * Bt ./ w;
z = Bs' * (W .* Z) * Bt ./ w;


end


function u = b(t,T,j,p)
i = j+1;
if p==0
    u = (t>= T(i) & t < T(i+p+1)); return 
end

u = zeros(size(t));
if T(i) ~= T(i+p)
    u = u + ((t-T(i)) / (T(i+p) -T(i))) .* b(t,T,j,p-1);
end
if T(i+1) ~= T(i+p+1)
    u = u + ((T(i+p+1)-t) ./ (T(i+p+1) -T(i+1))) .* b(t,T,j+1,p-1);
end
end





