function[fitresult, goodnessOfFit, xData, yData] = deconvCreateFit(xAxisConv, experimentalDataTemp, IRF, fitCoefficientsTwo)

%=====Preparing the Curve Data for Fitting====
%standardizing size, real numbers,removing NaN/Inf, nondoubles to doubles
[xData, yData, zData] = prepareCurveData(xAxisConv, experimentalDataTemp, IRF);

%======Performing fitting=====
fittingMethod = fittype('convolutionFitBiexponential(x,z,a,b,c,d)','dependent', {'y'}, 'independent', {'x', 'z'}, 'coefficients', {'a','b','c','d',}); 
fittingOptions = fitoptions( 'Method', 'NonlinearLeastSquares', 'MaxIter', 1000000000, 'MaxFunEval', 10000000, ...
    'StartPoint', fitCoefficientsTwo, 'TolX', 0.000000000001);
[fitresult, goodnessOfFit] = fit( [xData, zData], yData, fittingMethod, fittingOptions);

coefficients = coeffvalues(fitresult);
a = coefficients(1);
b = coefficients(2);
c = coefficients(3);
d = coefficients(4);
figure

tempLength = length(xData);
finalConv = conv(a*exp(b*xData)+c*exp(d*xData),IRF);
finalConv = finalConv(end-tempLength+1:end);
plot(xAxisConv,finalConv);
hold on
plot(xAxisConv, yData);
plot(xAxisConv, zData);

end

