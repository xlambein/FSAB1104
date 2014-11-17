function test_matlab4()

errReq  = 0.001;
Ih = gaussIntegrate(15,@z,errReq);
I = dblquad(@(x,y) z(x,y),-15,15,-15,15);
errObs = abs(I-Ih);

fprintf('\n\nComputation by Gauss-Legendre composite method \n');
fprintf('  Volume (Error required-observed) \n');
fprintf('  %5.7f (%1.7f-%1.7f)\n',Ih, errReq,errObs);

end

function z= z(x,y)
z = 8+19*exp(-0.005*(2*x.^2+y.^2));
end


