function [newData1] = importspectra(fileToRead1)
%IMPORTFILE1(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read
 
%  Auto-generated by MATLAB on 07-Jul-2022 11:12:22
 
DELIMITER = '\t';
HEADERLINES = 1;

% Import the file
newData1 = importdata(fileToRead1, DELIMITER, HEADERLINES);

