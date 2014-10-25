function test_matlab3()

n = 9; 
X = linspace(0,1,n+1);
U = exp(-2*X).*sin(10*pi*X);
dU = 10*pi*exp(-2*X).*cos(10*pi*X) -2 *U;

x = linspace(0,1,100);
u = exp(-2*x).*sin(10*pi*x);

tic
uh = hermite(n,X,U,dU,x);
t1=toc;

tic
hermiteslow(n,X,U,dU,x);
toc-t1

close all;
plot(x,u,'b-', x, uh, 'r-', X,U,'b.','Markersize',25);
fprintf('==== uh(18) %14.7e \n',uh(18));
 
end
