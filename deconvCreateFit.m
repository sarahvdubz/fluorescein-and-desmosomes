function[fitresult, goodnessOfFit, xData, yData] = deconvCreateFit(xAxisConv, experimentalDataTemp, IRF, fitCoefficients)

%=====Preparing the Curve Data for Fitting====
%standardizing size, real numbers,removing NaN/Inf, nondoubles to doubles
[xData, yData, zData] = prepareCurveData(xAxisConv, experimentalDataTemp, IRF);

%======Performing fitting=====
fittingMethod = fittype('convolutionFit(x,z,a,b)','dependent', {'y'}, 'independent', {'x', 'z'}, 'coefficients', {'a','b'}); 
fittingOptions = fitoptions( 'Method', 'NonlinearLeastSquares', 'MaxIter', 1000000000, 'MaxFunEval', 10000000, ...
    'StartPoint', fitCoefficients, 'TolX', 1e-19, 'TolFun', 1e-20);
[fitresult, goodnessOfFit] = fit( [xData, zData], yData, fittingMethod, fittingOptions);

coefficients = coeffvalues(fitresult);
a = coefficients(1);
b = coefficients(2);
figure

tempLength = length(xData);
finalConv = conv(a*exp(b*xData),IRF);
finalConv = finalConv(end-tempLength+1:end);
plot(xAxisConv,finalConv);
hold on
plot(xAxisConv, yData);
plot(xAxisConv, zData);

end

