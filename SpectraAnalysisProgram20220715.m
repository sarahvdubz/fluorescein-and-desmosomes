%================INSTRUCTIONS FOR USE================
%This program is used to plot multiple spectra on the same figure for 
%convenient comparison 

%=====Retrieving Files=====
% Specify the folder where the files live.
myFolder = 'C:\Users\Sarah Van Winkle\MATLAB Drive\20220715_fluorescein\Fluorescein Spectra';
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.dat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);



for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', fullFileName)
    data = importspectra(fullFileName);

%=====Normalizing Counts======
max_value = max(data.data(:,2));
min_value = min(data.data(:,2));
totalSpectra = data.data(:,2) - min_value;
maxTotalSpectra = max(totalSpectra);
totalSpectra = totalSpectra/maxTotalSpectra;




%===Generating X-Axis (assumes uniform spectra lengths)=====
x_axis = data.data(1:1024)


%===Generating Figure====
figure(2)
%add more colors as needed
colors = {'r','g','b','c','m','y','k'}
plot(x_axis, totalSpectra,'color',colors{k})
ylabel('Normalized Counts')
xlabel('Wavelength (nm)')
title('Emission Spectra')

hold on
end

labels = {"0% Glycerol", "20% Glycerol", "40% Glycerol", "60% Glycerol", "80% Glycerol"};
axis([450 700 0 1]);
legend(labels);