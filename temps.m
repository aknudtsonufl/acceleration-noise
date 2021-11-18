clc
clear all

filename = "./GraceFO Data/AHK1B_2020-01-01_D_04.txt";
temp = processAHK1B(filename);
disp(temp)

%% Data Loaded
% x = [];
% y = [];
tstruct = struct;

for ii=1:length(temp)
    tstruct(ii).x = temp(ii).x; %time in seconds
    tstruct(ii).y = temp(ii).y; %temp in C
end

%z+ = 3
%z- = 2
%y+ = 6
%y- = 4
% dy = tstruct(2).y-tstruct(3).y;
% [pxx1,f1] = pwelch(dy,hanning(length(dy)),[],[], 1/5);
% dT1 = sqrt(pxx1)'/(1e-3/f1).^(1/3);
t = temp(3)-temp(2);
%% psd
p = sqrt(psd(t));
t = sqrt(psd(temp(6)-temp(4)));

%% Plot
iplot(p,t)
iplot(sqrt(psd(temp([2;3;4;6]))))
% figure
% loglog(f1,sqrt(pxx1));
% 
% figure
% plot(tstruct(3).x, tstruct(3).y-tstruct(2).y)
% % xlabel('Time (s)')
% % ylabel('Temperature (C)')

%% Temp approx
% http://www.braeunig.us/space/atmmodel.htm
Tinf = 1000;
T10 = 360;
lamd = 0.01875;
z = 500;
z10 = 120;
ro = 6356.766;
eps = (z-z10)*(ro+z10)/(ro+z);



T = Tinf - (Tinf - T10)*exp(-lamd*eps)