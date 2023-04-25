%=================   Lifetime Analysis Tool v2     =====================

% DESCRIPTION: This program loads your raw TRPL .dat files, performs a
% deconvolution, and provides the final calculated lifetime.

% INSTRUCTIONS: Place this program within the folder containing your
% lifetime .dat files. It will load in ALL .dat files.

% ADDITIONAL PROCESSING:
%   1) Glycerol graph program: glycerolLifetimeGraph.m

% USER SETTINGS:
%Number of data channels?
dataChannelNumber = 2;

% Initializing Variables
allFileNames = [];
allFileNamesNoIRF = [];
allFilePaths = {};
allMaxIndex = [];
allMinIndex = [];
allNoiseIndex = [];
allFitCoeff = [];
allFitCoeff2 = [];
allFitCoeff3 = [];
allFitResult = [];
allDeconvGoF = [];
experimentalData = [];
goodnessOfFit = [];
goodnessOfFit2 = [];
goodnessOfFit3 = [];

% ====== STEP 1: Load experimental .dat files & extract TRPL data =======
currentDirectory = pwd;
fileSearchPattern = fullfile(currentDirectory, '*.dat'); 
fileData = dir(fileSearchPattern);

for k = 1:length(fileData)
    baseFileName = fileData(k).name;
    fullFilePath = fullfile(fileData(k).folder, baseFileName);
    allFilePaths{k} = fullFilePath;
    stringBaseFileName = string(baseFileName);
    allFileNames = [allFileNames; stringBaseFileName];
    fprintf(1, 'Now loading %s\n', fullFilePath);
end

%Finds IRF .dat files and isolates them
irfSearch = strfind(allFileNames, 'IRF');

for k = 1:length(allFilePaths)
    if isempty(irfSearch{k}) == 0;
        irfLocation = allFilePaths{k};
        irfIndex = k;
        allFilePaths{k} = [];
    else
        allFileNamesNoIRF = [allFileNamesNoIRF; allFileNames(k)];
    end
end

%allFileNames without IRF (for table summary)


allFilePaths = allFilePaths(~cellfun(@isempty,allFilePaths));
irfStruct = importfile(irfLocation);
IRF = irfStruct.data(:,1) + irfStruct.data(:,2);

for k = 1:length(allFilePaths)
    experimentalCellStruct{k} = importfile(allFilePaths{k});
    experimentalStruct = experimentalCellStruct{k};
    experimentalDataAll{k} = experimentalStruct.data(:,1) + experimentalStruct.data(:,2);
    experimentalHeaderAll{k} = experimentalStruct.textdata;
end

%Retrieving bin number and width
experimentalTextData = experimentalHeaderAll{1};
binSize = cell2mat(experimentalTextData(9));
binSize = extractBefore(binSize, 7);
binSize = str2double(binSize);
binNumber = cell2mat(experimentalTextData(3));
binNumber = str2double(binNumber);



%========= STEP 2: Trimming from Max to Noise Floor ======
%Cuts everything before maximum & normalizes data
for k = 1:length(experimentalDataAll)
    [maxLifetime, maxLifetimeIndex] = max(experimentalDataAll{k});
    allMaxIndex = [allMaxIndex; maxLifetimeIndex];
    normLifetime = experimentalDataAll{k}/maxLifetime;
    experimentalDataAllMax{k} = normLifetime(maxLifetimeIndex:end);
end

% Finds min positions 
for k = 1:length(experimentalDataAllMax)
    [minLifetime, minLifetimeIndex] = min(experimentalDataAllMax{k});
    allMinIndex = [allMinIndex; minLifetimeIndex];
end

experimentalDataLength = max(allMinIndex);

% trims to longest min position so they all match lengths
for k = 1:length(experimentalDataAllMax)
    normMaxLifetime = experimentalDataAllMax{k};
    experimentalDataAllMaxToMin{k} = normMaxLifetime(1:experimentalDataLength);
end

%========= STEP 3: Create Time Axis & IRF =======
xMax = experimentalDataLength * binSize;
xAxis = 0:binSize:xMax-binSize;
xAxis = xAxis';
[maxIRF, maxIRFIndex] = max(IRF);
convIRF = IRF;
convIRF = IRF(maxIRFIndex:maxIRFIndex+1651)/maxIRF;


%========= STEP 4: Preliminary Exponential Fitting for StartPoints =======
for k = 1:length(experimentalDataAllMaxToMin)

%===One-Term Exponential Fitting Tool===
[fitresult, gof2, xData, yData] = createFit(xAxis, ...
    experimentalDataAllMaxToMin{k});
fitCoefficients = coeffvalues(fitresult);
allFitCoeff = [allFitCoeff;fitCoefficients];
goodnessOfFitTemp = gof2.rsquare;
goodnessOfFit = [goodnessOfFit;goodnessOfFitTemp];

%===Two-Term Exponential Fitting Tool===
[fitresult, gof2, xData, yData] = createFitTwoTerms(xAxis, ...
    experimentalDataAllMaxToMin{k});
fitCoefficientsTwo = coeffvalues(fitresult);
allFitCoeff2 = [allFitCoeff2;fitCoefficientsTwo];
goodnessOfFitTemp2 = gof2.rsquare;
goodnessOfFit2 = [goodnessOfFit2;goodnessOfFitTemp2];

%===Three-Term Exponential Fitting Tool===
[fitresult, gof2, xData, yData] = createFitThreeTerms(xAxis, ...
    experimentalDataAllMaxToMin{k});
fitCoefficientsThree = coeffvalues(fitresult);
allFitCoeff3 = [allFitCoeff3;fitCoefficientsThree];
goodnessOfFitTemp3 = gof2.rsquare;
goodnessOfFit3 = [goodnessOfFit3;goodnessOfFitTemp3];

end

%=====STEP 4b: Calculating Preliminary Lifetimes=====
lifetimes = allFitCoeff(:,2);
lifetimes = (1./lifetimes) * -1;

biexponentialLifetimes = [allFitCoeff2(:,2), allFitCoeff2(:,4)]; 
biexponentialLifetimes = (1./biexponentialLifetimes) * -1;

triexponentialLifetimes = [allFitCoeff3(:,2), allFitCoeff3(:,4), allFitCoeff3(:,6)]; 
triexponentialLifetimes = (1./triexponentialLifetimes) * -1;

%====== STEP 5: Deconvolution with IRF ========

for k = 1:length(experimentalDataAllMaxToMin)
    inputFitCoeff = allFitCoeff(k,:);
    [tempFitResult, tempGoF] = deconvCreateFit(xAxis, experimentalDataAllMaxToMin{k},convIRF,  inputFitCoeff);
    tempFitResultCoeff = coeffvalues(tempFitResult);
    allFitResult = [allFitResult;tempFitResultCoeff];
    tempRSquare = tempGoF.rsquare;
    allDeconvGoF = [allDeconvGoF; tempRSquare];
end

%=====STEP 6: Final Deconvoluted Lifetime Calculation =====
lifetimesDeconv = allFitResult(:,2);
lifetimesDeconv = (1./lifetimesDeconv) * -1;
deconvolutedLifetimes = lifetimesDeconv(1:end);



%======STEP 7: Final Program Outputs for User====
%Command Window Output
fprintf(1, 'Displaying summary:');
summary = table(allFileNamesNoIRF, deconvolutedLifetimes, lifetimes, biexponentialLifetimes, triexponentialLifetimes)
save('lifetimeAnalysisVariables.mat');
