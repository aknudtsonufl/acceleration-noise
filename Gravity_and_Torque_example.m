clc
clear all

%% Setup
F = [];
T = [];
Grad = [];
dict = [];
m = []; %some nx1 mass array
r = []; %some nx3 position array

%% Calculations

%%%%%%%%%%Gravity from onboard masses%%%%%%%%%%
% F = [];
% T = [];
for ii=1:size(m)
    [grav_F,grav_T,grav_grad] = Gravity_ForceandTorque(m(ii),r(ii,:));
    F = [F; grav_F];
    T = [T; grav_T];
    Grad = [Grad; grav_grad];
    dict = [dict; strcat("point mass ",int2str(ii))];
end
disp('end')
% plot3(r(:,1),r(:,2),r(:,3),'*')
% grid()