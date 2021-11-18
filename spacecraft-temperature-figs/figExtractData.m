clc
clear all
uiopen('Laser-Temperature-telemetry-Part-A.fig', 1)
h = get(gca);
lines = h.Children;

x1 = lines(1).XData*86400;
y1 = lines(1).YData;
x2 = lines(2).XData*86400;
y2 = lines(2).YData;
x3 = lines(3).XData*86400;
y3 = lines(3).YData;
x4 = lines(4).XData*86400;
y4 = lines(4).YData;

[pxx1,f1] = pwelch(y1, [],[],[], 1/16);
[pxx2,f2] = pwelch(y2, [],[],[], 1/16);
[pxx3,f3] = pwelch(y3, [],[],[], 1/16);
[pxx4,f4] = pwelch(y4, [],[],[], 1/16);


figure
loglog(f1,sqrt(pxx1),f2,sqrt(pxx2),f3,sqrt(pxx3),f4,sqrt(pxx4))
grid()
xlabel('Frequency (Hz)')
ylabel('Temperature Spectrum (K/Hz^{1/2})')

dT1 = sqrt(pxx1)'/(1e-3/f1).^(1/3);
dT2 = sqrt(pxx2)'/(1e-3/f2).^(1/3);
dT3 = sqrt(pxx3)'/(1e-3/f3).^(1/3);
dT4 = sqrt(pxx4)'/(1e-3/f4).^(1/3);

%% ahk file
filename = "./GraceFO Data/AHK1B_2020-01-01_D_04.rpt";
temps = processAHK1B(filename);