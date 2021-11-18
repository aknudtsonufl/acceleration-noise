function [F,T] = Outgassing_ForceandTorque(dT)
    %dF_OGdT = A_wall*A*Io*theta/(T0^2)/Ceff; %N/K
    dF_OGdT = 5.7e-13; %N/K
    F = [dF_OGdT*dT,0,0];
    T = [0,0,0];
end