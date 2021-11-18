clc
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



%% Clear temporary variables

clearvars filename startRow formatSpec fileID dataArray ans;




%% Inputs
m = [1; 2]; r = [1,0,0; 0,1,0];          %Point masses and position relative to TM center, kg; m
s = 0.03; A = s^2; %m, m^2
T0 = 293.15; dT = 0.005; %K
P = 1e-5; %Pa
S = [3 2 0.8];

%positions
micro_rang = [-0.452, 0, -0.079];
uso1 = [-0.474, 0.358, 0];
uso2 = [-0.474, 0.358, 0];
ipu1 = [-0.387, 0.358, -0.079];
ipu2 = [-0.387, -0.358, -0.079];
gps_nav_ant = [-0.124523507	0	0.664556962];
r_micro = [micro_rang; uso1; uso2; ipu1; ipu2; gps_nav_ant];

lri_bat = [-0.318, 0.271, 0];
lri_notsure = [-0.145, 0.358, 0];
lri = [-0.343, 0, 0];
lri_baf = [-0.343, 0, 0];
lri_brac1 = [-0.110, 0.217, -0.079];
lri_brac2 = [-0.110, -0.217, -0.079];
lri_oth2 = [-0.118, -0.378, -0.079];
r_lri = [lri_bat; lri_notsure; lri; lri_baf; lri_brac1; lri_brac2; lri_oth2];

st1 = [0.0952, 0, 0];
tank1 = [-0.215, 0, 0];
tank2 = [0.198, 0, 0];
tank1m = [-0.215, 0, 0];
tank2m = [0.198, 0, 0];
r_attcont = [st1; tank1; tank2; tank1m; tank2m];

baf_end = [0.321, 0, -0.079];
black_box_end = [0.40343075	0	0.433544304];
black_box_upper = [0.390724269, 0.358041958, 0.21835443];
black_box_cl = [0.340533672	-0.123076923	0.21835443];
black_box_in = [0.244599746	0.267132867	0.158227848];
sol_bat_u = [0.045108005	0.281118881	0.205696203];
sol_bat_l = [0.045108005	-0.281118881	0.205696203];
SP1 = [0	0.366433566	0.458860759];
SPt = [0	0	0.683544304];
SP2 = [0	-0.366433566	0.458860759];
r_misc = [baf_end; black_box_end; black_box_upper; black_box_cl; black_box_in; sol_bat_u; sol_bat_l; SP1; SPt; SP2];

r = [r_micro
    r_lri
    r_attcont
    r_misc];
r = r.*[3 2 0.8];

%volumes
micro_rangV = 0.0104;
uso1V = 0.002023;
uso2V = 0.002023;
ipu1V = 0.002231;
ipu2V = 0.002231;
gps_nav_antV = 0.001279531;
V_micro = [micro_rangV; uso1V; uso2V; ipu1V; ipu2V; gps_nav_antV];

lri_batV = 0.001254;
lri_notsureV = 0.002295;
lriV = 0.005774/2;
lri_bafV = 0.005774/2;
lri_brac1V = 0.000514;
lri_brac2V = 0.000514;
lri_oth2V = 0.001939;
V_lri = [lri_batV; lri_notsureV; lriV; lri_bafV; lri_brac1V; lri_brac2V; lri_oth2V];

st1V = 0.003689;
tank1V = 0.006795;
tank2V = 0.006795;
tank1m = 0.006795;
tank2m = 0.006795;
V_attandcont = [st1V; tank1V; tank2V; tank1m; tank2m];

baf_endV = 0.005774/2;
black_box_endV = 0.014969;
black_box_upperV = 0.00178504;
black_box_clV = 0.00214946;
black_box_inV = 0.003361643;
sol_bat_uV = 0.007640558;
sol_bat_lV = 0.007640558;
SP1V = 0.073628408;
SPtV = 0.096220968;
SP2V = 0.073628408;
%include the electrode housing - makes a difference in the gradient
V_misc = [baf_endV; black_box_endV; black_box_upperV; black_box_clV; black_box_inV; sol_bat_uV; sol_bat_lV; SP1V; SPtV; SP2V];


V = S(1)*S(2)*S(3)*[V_micro;
    V_lri;
    V_attandcont;
    V_misc];
V(17:18) = V(17:18)-4/3*pi*((V(17:18)*3/4/pi).^(1/3)-0.005).^3;

m = V*1200;
m(7) = 5.1;
m(8) = 15.5;
m(9) = 1.1/3;
m(10:11) = m(9:10);
m(14:15) = m(13:14)/1200*480;%300;
m(17:18) = m(17:18)/1200*4500;
m(19) = 1.1/3;
m_end = length(m);
m(m_end-2:m_end) = m(m_end-2:m_end)/0.03/1200*2.06;
pos_r = [];
neg_r = [];
for ii=1:length(r)
    if(r(ii)>0)
        pos_r = [pos_r, ii];
    elseif(r(ii)<0)
        neg_r = [neg_r, ii];
    end
end

pos_test = m(pos_r)./r(pos_r,1).^2;
neg_test = m(neg_r)./r(neg_r,1).^2;
%Calculate gradient of gravity force

%% Constants
sigma = 5.67e-8;                %stephan boltzmann constant (W*m^-2*K^-4)
c  = 3e8;                       %speed of light (m/s)
k_B = 1.38064852e-23;           %Boltzmann constant (m2 kg s-2 K-1)


%% Setup
F = [];
T = [];
Grad = [];
dict = [];

%% Calculations

%%%%%%%%%%Gravity from onboard masses%%%%%%%%%%
for ii=1:size(m)
    [grav_F,grav_T,grav_grad] = Gravity_ForceandTorque(m(ii),r(ii,:));
    F = [F; grav_F];
    T = [T; grav_T];
    Grad = [Grad; grav_grad];
    dict = [dict; strcat("point mass ",int2str(ii))];
end

% plot3(r(:,1),r(:,2),r(:,3),'*')
% grid()
%% Magnetic onboard
X0 = 2e-5;
mu0 = pi*4e-7; %magnetic field constant
sig0 = 1;
Nd = 1;
Nf = 1000;
f = logspace(-4, 0, Nf)';
f = f(:, ones(1, Nd));
X_tm = 2e-5;%X0-j*s*mu0*sig0*pi*f/12;
M = 1e-9; %check this number

dipoles = zeros(1,6); %First 3 columns give x,y,z location of dipole, 4-6th gives magnetic moment

dipoles = [-0.2570 -0.6247 0.5904 261.969 261.969 261.969   %Caging Mechanicsm Control Unit
            -0.4238 -0.6285 0.5989 38.033 38.033 38.033     %Data Management Unit Dipole 1
            -0.3814 -0.5029 0.5750 38.033 38.033 38.033     %Data Management Unit Dipole 2
            -0.2506 -0.3325 0.4139 0.200 0.200 0.200        %IS FEE PCU Dipole 1
            -0.3838 -0.5211 0.2869 19.119 19.119 19.119     %IS FEE PCU Dipole 2
            -0.4111 -0.5804 0.3019 6.235 6.235 6.235        %IS FEE PCU Dipole 3
            -0.3021 -0.6456 0.2539 7.063 7.063 7.063        %IS FEE PCU Dipole 4
            -0.5037 0.6208 0.6973 5.633 5.633 5.633         %IS FEE SAU Dipole 1
            -0.4579 0.4856 0.4502 12.001 12.001 12.001      %IS FEE SAU Dipole 2
            -0.2675 0.7150 0.4536 5.633 5.633 5.633         %IS FEE SAU 2 Dipole 1
            -0.1732 0.6076 0.7005 12.001 12.001 12.001      %IS FEE SAU 2 Dipole 2
            0.7873 -0.7460 0.0618 6.512 6.512 6.512         %Laser Assembly Dipole 1
            0.5140 -0.7826 0.1578 8.958 8.958 8.958         %Laser Assembly Dipole 2
            0.6049 -0.6400 0.2089 17.295 17.295 17.295      %Laser Assembly Dipole 3
            0.5659 -0.6025 0.0169 1.459 1.459 1.459         %Laser Assembly Dipole 4
            0.6486 -0.6908 0.5478 3.362 3.362 3.362         %Laser Modulation Unit Dipole 1
            0.6777 -0.6984 0.5218 3.977 3.977 3.977         %Laser Modulation Unit Dipole 2
            0.5660 -0.2188 0.5488 77.455 77.455 77.455      %Phasemeter Unit Dipole 1
            0.5890 -0.2998 0.6338 60.814 60.814 60.814      %Phasemeter Unit Dipole 2
            -0.5517 -0.1842 0.3591 20.151 20.151 20.151     %Radiation Monitor
            0.6366 -0.7295 0.2679 4.236 4.236 4.236         %Reference Laser Unit Dipole 1
            0.5960 -0.7419 0.3069 6.129 6.129 6.129         %Reference Laser Unit Dipole 2
            -0.3887 -0.8187 1.0107 1.158 1.158 1.158        %Sun Sensor 1
            0.3906 0.9740 1.0507 0.592 0.592 0.592          %Sun Sensor 2
            0.5106 0.8577 1.0107 1.225 1.225 1.225          %Sun Sensor 3
            -0.4258 -0.3368 0.5475 35.157 35.157 35.157     %Ultra-Violet Lamp Unit
            -0.2170 -0.8650 0.5560 824.621 824.621 824.621  %Dipole 1
            -0.7620 0.4120 0.5200 824.621 824.621 824.621   %Dipole 2
             0.6330 0.6240 0.5590 824.621 824.621 824.621   %Dipole 3      
    ];
un = random_unit_vector(length(dipoles))';
dipoles(:,4:6) = dipoles(:,4:6).*un;
dipoles(:,1:3) = dipoles(:,1:3).*S;


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
clear MAG1B20180610D04
clear GF1AccTemps
clear GFAccTemps

Be = arr(:,3)/1e6;
Bd = zeros(length(Be),3);
Bgrad = zeros(length(Be),3);
Bgrad_d = zeros(length(Be),3);
Bgrad_e = zeros(length(Be),1);
x = [0,0,0];

for a=1:length(Be) %give time series
    for b=1:length(dipoles(:,1)) %Question: information given is magnitude of magnetic moment, equation is looking for vector?
        xa = dipoles(b,1:3);
        ma = dipoles(b,4:6);
        rhat = (x-xa)/norm(x-xa);
        Bda = mu0/4/pi*(3*dot(ma,rhat)*rhat-ma)/norm(x-xa)^3; % Equation from page 91 of https://upcommons.upc.edu/bitstream/handle/2117/95156/TMDA1de1.pdf?sequence=1&isAllowed=y
        Bd(a,:) = Bd(a,:)+Bda;
        Bgrad_d_a = -3*Bd(a,:)/norm(x-xa);
        Bgrad_d(a,:) = Bgrad_d(a,:)+Bgrad_d_a;
    end
    Bgrad_e(a) = -3*Be(a)/(495000+earthRadius);
end

%% Force calc
B = Be+Bd(:,1);
Bgrad = Bgrad_d(:,1);%+Brad_e;
% Fbe = (X_tm*s^3/mu0.*Be+M).*Bgrad_e;
% Fbd = (X_tm*s^3/mu0.*Bd(:,1)+M).*Bgrad_d(:,1);
Fbx = (X_tm*s^3/mu0.*B+M).*Bgrad;
%% Plot
figure
plot(t,Fbx(:,1),t,Fbe,'--',t,Fbd,'-.')
grid()
xlabel('Time')
ylabel('Force (N)')
legend('Total','Earth','Dipoles')

figure
plot(t,Be,t,Bd(:,1),'--')
grid()
title('Magnetic Field')
xlabel('Time')
ylabel('microTesla')
legend('Earth Magnetic Field','Onboard Magnetic Field')

figure
plot(t,Bgrad_e,t,Bgrad_d(1),'--')
grid()
title('Magnetic Field Gradient')
xlabel('Time')
ylabel('microTesla/m')
legend('Earth Magnetic Field','Onboard Magnetic Field')

%%%%%%%%%%Accounting for Thermal Expansion%%%%%%%%%%
% F = zeros(length(GFAccTemps)-1,length(m));
% F_expsum = zeros(length(GFAccTemps)-1,1);
% for ii=2:length(GFAccTemps)
%     %disp(2700-ii)
%     dT = GFAccTemps(ii,10)-GFAccTemps(ii-1,10);
%     alpha = 1e-6;%(-1.5e-6*(T-23)^3+3.17e-4*(T-23)^2-2.44e-2*(T-23)+0.015)*10^-6;
%     r_exp = r+r.*alpha*dT;
%     F_exp = zeros(length(m),3);
%     T_exp = zeros(length(m),3);
%     for jj=1:size(m)
%         [grav_F,grav_T] = Gravity_ForceandTorque(m(jj),r_exp(jj,:));
%         F_exp(jj,:) = grav_F;
%         T_exp(jj,:) = grav_T;
%         dict = [dict; strcat("point mass ",int2str(jj))];
%     end
%     F(ii-1,:) = F_exp(:,1);
%     F_expsum(ii-1) = sum(F_exp(1));
% end
% F_diff = F_expsum(2:2699)-F_expsum(1:2698);
% 
% 
% %% Run separately
% [pxx1,f1] = pwelch(F_diff, [],[],[], 1/32);
% figure(2)
% loglog(f1,sqrt(pxx1))
% xlabel('Frequency (Hz)')
% ylabel('N/Hz^{1/2}')
% grid()

%%%%%%%%%%Thermal Effects%%%%%%%%%%

% %Radiometer
% [radio_F, radio_T] = Radiometer_ForceandTorque(A,P,T0,dT);
% F = [F; radio_F];
% T = [T; radio_T];
% dict = [dict; "Radiometer Thermal"];
% 
% 
% %Radiation Pressure
% [rad_F, rad_T] = RadiationPressure_ForceandTorque(A,T0,sigma,c,dT);
% F = [F; rad_F];
% T = [T; rad_T];
% dict = [dict; "Radiation Pressure Thermal"];
% 
% 
% %Outgassing
% [outg_F, outg_T] = Outgassing_ForceandTorque(dT);
% F = [F; outg_F];
% T = [T; outg_T];
% dict = [dict; 'Outgassing Thermal'];
% 
%%%%%%%%%% %%%%%%%%%%

%%%%%%%%%% %%%%%%%%%%

%%%%%%%%%% %%%%%%%%%%

%%%%%%%%%%Total%%%%%%%%%%
% 
% %positions
% micro_rang = [-0.452, 0, -0.079];
% uso1 = [-0.474, 0.358, 0];
% uso2 = [-0.474, 0.358, 0];
% ipu1 = [-0.387, 0.358, -0.079];
% ipu2 = [-0.387, -0.358, -0.079];
% r_micro = [micro_rang; uso1; uso2; ipu1; ipu2];
% 
% lri_bat = [-0.318, 0.271, 0];
% lri_notsure = [-0.145, 0.358, 0];
% lri = [-0.343, 0, 0];
% lri_baf = [-0.343, 0, 0];
% lri_brac1 = [-0.110, 0.217, -0.079];
% lri_brac2 = [-0.110, -0.217, -0.079];
% lri_oth2 = [-0.118, -0.378, -0.079];
% r_lri = [lri_bat; lri_notsure; lri; lri_baf; lri_brac1; lri_brac2; lri_oth2];
% 
% st1 = [0.0952, 0, 0];
% tank1 = [-0.215, 0, 0];
% tank2 = [0.198, 0, 0];
% tank1m = [-0.215, 0, 0];
% tank2m = [0.198, 0, 0];
% r_attcont = [st1; tank1; tank2; tank1m; tank2m];
% 
% baf_end = [0.321, 0, -0.079];
% black_box = [0.045, -0.267, 0];
% black_box_end = [0.40343075	0	0.433544304];
% SP1 = [0	0.366433566	0.458860759];
% SPt = [0	0	0.683544304];
% SP2 = [0	-0.366433566	0.458860759];
% 
% r = [r_micro
%     r_lri
%     r_attcont
%     baf_end
%     black_box
%     black_box_end
%     SP1
%     SPt
%     SP2];
% r = r.*[3 2 0.8];
% 
% %volumes
% micro_rangV = 0.0104;
% uso1V = 0.002023;
% uso2V = 0.002023;
% ipu1V = 0.002231;
% ipu2V = 0.002231;
% V_micro = [micro_rangV; uso1V; uso2V; ipu1V; ipu2V];
% 
% lri_batV = 0.001254;
% lri_notsureV = 0.002295;
% lriV = 0.005774/2;
% lri_bafV = 0.005774/2;
% lri_brac1V = 0.000514;
% lri_brac2V = 0.000514;
% lri_oth2V = 0.001939;
% V_lri = [lri_batV; lri_notsureV; lriV; lri_bafV; lri_brac1V; lri_brac2V; lri_oth2V];
% 
% st1V = 0.003689;
% tank1V = 0.006795;
% tank2V = 0.006795;
% tank1m = 0.006795;
% tank2m = 0.006795;
% V_attandcont = [st1V; tank1V; tank2V; tank1m; tank2m];
% 
% baf_endV = 0.005774/2;
% black_boxV = 0.006703;
% black_box_endV = 0.014969;
% SP1V = 0.073628408;
% SPtV = 0.096220968;
% SP2V = 0.073628408;
% 
% 
% V = 3*2*0.8*[V_micro
%     V_lri
%     V_attandcont
%     baf_endV
%     black_boxV
%     black_box_endV
%     SP1V
%     SPtV
%     SP2V];
% V(16:17) = V(16:17)-4/3*pi*((V(15:16)*3/4/pi).^(1/3)-0.005).^3;
% 
% m = V*1200;
% m(7) = 5.1;
% m(8) = 15.5;
% m(9) = 1.1/3;
% m(10:11) = m(9:10);
% m(14:15) = m(13:14)/1200*480;%300;
% m(16:17) = m(15:16)/1200*4500;
% m(18) = 1.1/3;
% m(21:23) = m(21:23)/0.03/1200*2.06;
% 
