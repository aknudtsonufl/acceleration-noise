

%% Initialize variables.
% Specify the folder where the files live.
myFolder = 'C:\Users\aaronknudtson\Documents\ACC_TEMPERATURES\ACC_TEMPERATURES';
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.txt'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
GFAccTemps = [];
GFAccTimes = [];
for k = 1
    baseFileName = theFiles(k).name;
    filename = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', filename);
    
    startRow = 2;
    formatSpec = '%10{yyyy-MM-dd}D%9{HH:mm:ss}D%14f%16f%20f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%16f%f%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'EmptyValue', NaN, 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    fclose(fileID);
    GF1AccTemps = table(dataArray{1:end-1}, 'VariableNames', {'YYYYMMDD','HHMMDD','DOY','MJD','tSDS','THT10042_enTH','T10052_enTHT10','enTHT10116','enTHT10043_en','THT10075_enT','HT10053_enTHT1','enTHT1011','enTHT10030_e','nTHT10025_en','THT10052_enTHT','enTHT101','enTHT10053_','enTHT10088_en','VarName19','VarName20','VarName21'});
    GFAccTemps = [GFAccTemps, table2array(GF1AccTemps(:,3:end))];
    GFAccTimes = [GFAccTimes, table2array(GF1AccTemps(:,2))];
end
count = [];
for ii=4:length(GFAccTemps(1,:))
    plot(GFAccTimes,GFAccTemps(:,ii))
    hold on
end
    
    




%% Clear temporary variables

clearvars filename startRow formatSpec fileID dataArray ans;