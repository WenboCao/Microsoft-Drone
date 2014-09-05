%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Users/orangemr/Documents/MATLAB/matlab/test.CSV
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2014/08/15 14:48:59

%% Initialize variables.
filename = '/Users/orangemr/Documents/MATLAB/matlab/test.CSV';
delimiter = ',';

%% Read columns of data as strings:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric strings to numbers.
% Replace non-numeric strings with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,3,4,6,8,9,10,11,12,13,14,16,17,18,19,20]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(thousandsRegExp, ',', 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end

%% Split data into numeric and cell columns.
rawNumericColumns = raw(:, [2,3,4,6,8,9,10,11,12,13,14,16,17,18,19,20]);
rawCellColumns = raw(:, [1,5,7,15,21]);


%% Allocate imported array to column variable names
VarName1 = rawCellColumns(:, 1);
VarName2 = cell2mat(rawNumericColumns(:, 1));
VarName3 = cell2mat(rawNumericColumns(:, 2));
VarName4 = cell2mat(rawNumericColumns(:, 3));
N = rawCellColumns(:, 2);
VarName6 = cell2mat(rawNumericColumns(:, 4));
W = rawCellColumns(:, 3);
VarName8 = cell2mat(rawNumericColumns(:, 5));
G3 = cell2mat(rawNumericColumns(:, 6));
VarName10 = cell2mat(rawNumericColumns(:, 7));
VarName11 = cell2mat(rawNumericColumns(:, 8));
VarName12 = cell2mat(rawNumericColumns(:, 9));
VarName13 = cell2mat(rawNumericColumns(:, 10));
VarName14 = cell2mat(rawNumericColumns(:, 11));
VarName15 = rawCellColumns(:, 4);
VarName16 = cell2mat(rawNumericColumns(:, 12));
VarName17 = cell2mat(rawNumericColumns(:, 13));
VarName18 = cell2mat(rawNumericColumns(:, 14));
VarName19 = cell2mat(rawNumericColumns(:, 15));
VarName20 = cell2mat(rawNumericColumns(:, 16));
VarName21 = rawCellColumns(:, 5);

%% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me rawNumericColumns rawCellColumns;