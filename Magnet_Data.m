%% Import data from text file
% Script for importing data from the following text file:
%
%    filename: C:\Users\aaronknudtson\Downloads\GraceFO1009\MAG1B_2018-06-10_D_04.txt
%
% Auto-generated by MATLAB on 14-Oct-2020 16:47:16

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 43);

% Specify range and delimiter
opts.DataLines = [92001, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["header", "VarName2", "Var3", "Var4", "VarName5", "VarName6", "VarName7", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var30", "Var31", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38", "Var39", "Var40", "Var41", "Var42", "Var43"];
opts.SelectedVariableNames = ["header", "VarName2", "VarName5", "VarName6", "VarName7"];
opts.VariableTypes = ["double", "double", "string", "string", "double", "double", "double", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["Var3", "Var4", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var30", "Var31", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38", "Var39", "Var40", "Var41", "Var42", "Var43"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var3", "Var4", "Var8", "Var9", "Var10", "Var11", "Var12", "Var13", "Var14", "Var15", "Var16", "Var17", "Var18", "Var19", "Var20", "Var21", "Var22", "Var23", "Var24", "Var25", "Var26", "Var27", "Var28", "Var29", "Var30", "Var31", "Var32", "Var33", "Var34", "Var35", "Var36", "Var37", "Var38", "Var39", "Var40", "Var41", "Var42", "Var43"], "EmptyFieldRule", "auto");

% Import the data
MAG1B20180610D04 = readtable("C:\Users\aaronknudtson\Downloads\GraceFO1009\MAG1B_2018-06-10_D_04.txt", opts);


%% Clear temporary variables
clear opts

arr = table2array(MAG1B20180610D04);
timearr = seconds(arr(:,1)+arr(:,2)/1e9);
t1 = datetime(2000,1,1);
t = t1+timearr;
mags = arr(:,3);
plot(t,mags)
xlabel('Time')
ylabel('microTesla')

