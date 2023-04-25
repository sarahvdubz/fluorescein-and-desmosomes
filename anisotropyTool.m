%=================   Anisotropy Tool     =====================

% DESCRIPTION: This program loads your raw TRPL .dat files, performs a
% deconvolution, and uses the anisotropy equation to build the anisotropy curve.

% INSTRUCTIONS: Place this program within the folder containing your
% lifetime .dat files. It will load in ALL .dat files.

%Don't normalize, add the g factors

%=====STEP 1: Importing Deconvoluted Intensities=====
load("lifetimeAnalysisVariables.mat");

paraIn1 = experimentalDataAllMaxToMin{2};
perpIn1 = experimentalDataAllMaxToMin{1};

anisotropy = (paraIn1-perpIn1)./(paraIn1+2*perpIn1);

figure(1)
semilogy(xAxis,anisotropy);

[fitresult, gof2, xData, yData] = createFit(xAxis, ...
    anisotropy);
fitCoefficients = coeffvalues(fitresult);
goodnessOfFitTemp = gof2.rsquare;
rotCorrel = 1/fitCoefficients(1,2)


paraIn2 = experimentalDataAllMaxToMin{5};
perpIn2 = experimentalDataAllMaxToMin{4};

anisotropy2 = (paraIn2-perpIn2)./(paraIn2+2*perpIn2);
figure(2)
plot(xAxis,anisotropy2)