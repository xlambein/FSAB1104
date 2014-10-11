function test_matlab1()

%
% -1- Simple test d'exécution !
%

close all;
n = 20;
m = 1000;
X = linspace(0,5,n+1);
U = 1 ./ (X + 1);

x = linspace(eps,5*(1-eps),m);
uh = lagrange_super(x,X,U);

figure;
plot(x,uh,X,U,'.','Markersize',20);
title('The super Lagrange interpolation');



%
% -2- Test d'efficacite calculatoire
%
%     Que se passe-t-il lorsque n > 60 ?
%     Pourquoi ?
%

fprintf('===== Efficiency tests with m = %3d :-)\n',m);
nSet = [5 10 20 40 60];
m = 1000; i = 0;
for n = nSet
    i = i+1;
    X = linspace(0,10,n+1);
    fprintf('  === Order of interpolation is %3d == \n',n);

    x = linspace(eps,10*(1-eps),m);
    U = 1 ./ (X + 1);
    fprintf('   == Naive Lagrange  :',n);    
    tic; uh = lagrange_naive(x,X,U); na(i) = toc; toc
    fprintf('   == Clever Lagrange :',n);
    tic; uh = lagrange(x,X,U); la(i) = toc; toc  
    fprintf('   == Super Lagrange  :',n);    
    tic; uh = lagrange_super(x,X,U); su(i) = toc; toc
 
    
end

%
% -3- Plots de l'efficacité du calcul
%     Pourquoi effectuer un plot du temps de calcul en log-log ?    
%

figure
loglog(nSet,na,'.-r',nSet,la,'.-b',nSet,su,'.-g');
xlabel('Order : n');
ylabel('CPU time');
legend('naive','clever','super');
title('Computational efficiency');

end