clc
clear all

%% Setup the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 6);

% Specify range and delimiter
opts.DataLines = [10, Inf];
opts.Delimiter = " ";

% Specify column names and types
opts.VariableNames = ["VarName1", "Var2", "Var3", "Var4", "VarName5", "VarName6"];
opts.SelectedVariableNames = ["VarName1", "VarName5", "VarName6"];
opts.VariableTypes = ["double", "string", "string", "string", "double", "string"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";
opts.ConsecutiveDelimitersRule = "join";
opts.LeadingDelimitersRule = "ignore";

% Specify variable properties
opts = setvaropts(opts, ["Var2", "Var3", "Var4", "VarName6"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var2", "Var3", "Var4", "VarName6"], "EmptyFieldRule", "auto");

% Import the data
GRACEnoshield = readtable("C:\Users\aaronknudtson\Documents\Research\MATLAB\Noise Model\COMSOL simulations\GRACE_noshield.txt", opts);
GRACEshieldopenthickness1 = readtable("C:\Users\aaronknudtson\Documents\Research\MATLAB\Noise Model\COMSOL simulations\GRACE_shield_open_thickness_1.txt", opts);
GRACEshieldopenthickness2 = readtable("C:\Users\aaronknudtson\Documents\Research\MATLAB\Noise Model\COMSOL simulations\GRACE_shield_open_thickness_2.txt", opts);
GRACEshieldopenthickness4 = readtable("C:\Users\aaronknudtson\Documents\Research\MATLAB\Noise Model\COMSOL simulations\GRACE_shield_open_thickness_4.txt", opts);
GRACEshieldopenthickness5 = readtable("C:\Users\aaronknudtson\Documents\Research\MATLAB\Noise Model\COMSOL simulations\GRACE_shield_open_thickness_5.txt", opts);

GRACEshieldclosedhole5thickness1 = readtable("C:\Users\aaronknudtson\Documents\Research\MATLAB\Noise Model\COMSOL simulations\GRACE_shield_closedhole5_thickness_1.txt", opts);
GRACEshieldclosedhole5thickness2 = readtable("C:\Users\aaronknudtson\Documents\Research\MATLAB\Noise Model\COMSOL simulations\GRACE_shield_closedhole5_thickness_2.txt", opts);

GRACEstlshield2 = readtable("C:\Users\aaronknudtson\Documents\Research\MATLAB\Noise Model\COMSOL simulations\GRACE_stlshield_2.txt", opts);
GRACEstpshield = readtable("C:\Users\aaronknudtson\Documents\Research\MATLAB\Noise Model\COMSOL simulations\GRACE_shieldSTP.txt", opts);


% Table2array
GRACEnoshieldarr = table2array(GRACEnoshield);
GRACEnoshieldarr = str2double(GRACEnoshieldarr);

GRACEshieldopenthickness1arr = table2array(GRACEshieldopenthickness1);
GRACEshieldopenthickness1arr = str2double(GRACEshieldopenthickness1arr);

GRACEshieldopenthickness2arr = table2array(GRACEshieldopenthickness2);
GRACEshieldopenthickness2arr = str2double(GRACEshieldopenthickness2arr);

GRACEshieldopenthickness4arr = table2array(GRACEshieldopenthickness4);
GRACEshieldopenthickness4arr = str2double(GRACEshieldopenthickness4arr);

GRACEshieldopenthickness5arr = table2array(GRACEshieldopenthickness5);
GRACEshieldopenthickness5arr = str2double(GRACEshieldopenthickness5arr);



GRACEshieldclosedhole5thickness1arr = table2array(GRACEshieldclosedhole5thickness1);
GRACEshieldclosedhole5thickness1arr = str2double(GRACEshieldclosedhole5thickness1arr);

GRACEshieldclosedhole5thickness2arr = table2array(GRACEshieldclosedhole5thickness2);
GRACEshieldclosedhole5thickness2arr = str2double(GRACEshieldclosedhole5thickness2arr);



GRACEstlshield2arr = table2array(GRACEstlshield2);
GRACEstlshield2arr = str2double(GRACEstlshield2arr);

GRACEstpshieldarr = table2array(GRACEstpshield);
GRACEstpshieldarr = str2double(GRACEstpshieldarr);



%% Clear temporary variables
clear opts

figure
plot(GRACEnoshieldarr(:,1), GRACEnoshieldarr(:,2), '.')
hold on
plot(GRACEshieldopenthickness1arr(:,1), GRACEshieldopenthickness1arr(:,2), '.')
hold on
plot(GRACEshieldopenthickness2arr(:,1), GRACEshieldopenthickness2arr(:,2), '.')
hold on
plot(GRACEshieldopenthickness4arr(:,1), GRACEshieldopenthickness4arr(:,2), '.')
hold on
plot(GRACEshieldopenthickness5arr(:,1), GRACEshieldopenthickness5arr(:,2), '.')
legend('no shield', '1 mm thickness','2 mm thickness', '4 mm thickness', '5 mm thickness')
xlabel('x [m]')
ylabel('log10(Magnetic Field Strength [T] )')
title('Open top and bottom, variable thickness')


maxLength = max([length(GRACEnoshieldarr), length(GRACEshieldopenthickness1arr)]);
xFit = 1:maxLength;
interpnoThickness = interp1(1:length(GRACEnoshieldarr), GRACEnoshieldarr, xFit);
interpThickness1 = interp1(1:length(GRACEshieldopenthickness1arr), GRACEshieldopenthickness1arr, xFit);
% 
% figure
% plot(interpnoThickness(:,1), interpnoThickness(:,2)-interpThickness1(:,2), '.')
% legend('shield vs. no shield difference')

%% f
figure
plot(GRACEnoshieldarr(:,1), GRACEnoshieldarr(:,2), '.')
hold on
plot(GRACEshieldopenthickness1arr(:,1), GRACEshieldopenthickness1arr(:,2), '.')
hold on
plot(GRACEshieldclosedhole5thickness1arr(:,1), GRACEshieldclosedhole5thickness1arr(:,2), '.')
xlabel('x [m]')
plot(GRACEshieldclosedhole5thickness2arr(:,1), GRACEshieldclosedhole5thickness2arr(:,2), '.')
xlabel('x [m]')

ylabel('log10(Magnetic Field Strength [T] )')
title('Hole size 5 mm, variable thickness')
legend('no shield','1 mm thick open shield', '1 mm thick shield w/ hole')


%% stl shield
figure
plot(GRACEnoshieldarr(:,1), GRACEnoshieldarr(:,2), '.')
hold on
plot(GRACEshieldopenthickness1arr(:,1), GRACEshieldopenthickness1arr(:,2), '.')
hold on
plot(GRACEshieldclosedhole5thickness1arr(:,1), GRACEshieldclosedhole5thickness1arr(:,2), '.')
hold on
plot(GRACEstlshield2arr(:,1), GRACEstlshield2arr(:,2), '.')
hold on
plot(GRACEstpshieldarr(:,1), GRACEstpshieldarr(:,2), '.')
xlabel('x [m]')
ylabel('Magnetic Field (T)')


ylabel('log10(Magnetic Field Strength [T] )')
title('Hole size 5 mm, variable thickness')
legend('no shield','1 mm thick open shield', '1 mm thick shield w/ 2.5 mm radius hole', 'STP shield', 'STP shield working')

