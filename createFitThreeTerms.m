function[fitresult, goodnessOfFit, xData, yData] = createFitThreeTerms(xAxisLifetimeForFit, normalizedLifetimeMaxToNoiseForFit)

%=====Preparing the Curve Data for Fitting====
%standardizing size, real numbers,removing NaN/Inf, nondoubles to doubles
[xData, yData] = prepareCurveData(xAxisLifetimeForFit, normalizedLifetimeMaxToNoiseForFit);


%======Performing fitting=====
fittingMethod = fittype('a*exp(b*x)+c*exp(d*x)+e*exp(f*x)','dependent', {'y'}, 'independent', {'x'}, 'coefficients', {'a','b','c','d','e','f'}); 
fittingOptions = fitoptions( 'Method', 'NonlinearLeastSquares', 'MaxIter', 10000000, 'Upper', [Inf, 0, Inf, 0, Inf, 0], 'StartPoint', [0.5 -0.5 0.5 -0.5 0.5 -0.5] );
[fitresult, goodnessOfFit] = fit( xData, yData, fittingMethod, fittingOptions);

