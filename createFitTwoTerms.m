function[fitresult, goodnessOfFit, xData, yData] = createFitTwoTerms(xAxisLifetimeForFit, normalizedLifetimeMaxToNoiseForFit)

%=====Preparing the Curve Data for Fitting====
%standardizing size, real numbers,removing NaN/Inf, nondoubles to doubles
[xData, yData] = prepareCurveData(xAxisLifetimeForFit, normalizedLifetimeMaxToNoiseForFit);


%======Performing fitting=====
%(two-term exponential is built into MATLAB, so opting to use existing 
%fittype.

fittingMethod = fittype('exp2');
fittingOptions = fitoptions( 'Method', 'NonlinearLeastSquares', 'MaxIter', 1000000  );
[fitresult, goodnessOfFit] = fit( xData, yData, fittingMethod, fittingOptions);



