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



%========= STEP 4: Trimming from Max to Noise Floor ======
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

%========= STEP 3: Create Time Axis =======
xMax = experimentalDataLength * binSize;
xAxis = 0:binSize:xMax-binSize;
xAxis = xAxis';

%===== STEP 5: Plot all curves ======
for k = 1:length(experimentalDataAllMaxToMin)
    plot(xAxis, experimentalDataAllMaxToMin{k})
    hold on
end
hold on
legend(allFileNamesNoIRF);