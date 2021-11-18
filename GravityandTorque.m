G = 6.67408e-11;                %Gravitational Constant
n = 10;                         %Number of Divisions
s = 0.03;                       %side length
m_TM = 0.54;                    %mass of TM
m = 1; r = [0, 0, 0];           %point mass and location
I   = 1/6*m*s^2;                %Moment of Inertia



F=0;
T=0;
rr = [0,0,0];
for ii=1:n
    for jj=1:n
        for kk=1:n
            r_c = [-s/2 + s*(2*ii-1)/(2*n), -s/2 + s*(2*jj-1)/(2*n), -s/2 + s*(2*kk-1)/(2*n)]; %center of mass of piece to center of mass of TM
            r_p = r+r_c; %Takes care of any offsets
%             rr = rr+r_p;
            F_i = G*m*m_TM/n^3/norm(r_p)^2.*r_p/norm(r_p);                  %second part is the directional vector
            T_i = cross(r_c,F_i);
            F = F + F_i;
            T = T + T_i;
        end
    end
end

% points = s/4*[1, 1, 1
%             1, 1, -1
%             1, -1, 1
%             1, -1, -1
%             -1, 1, 1
%             -1, 1, -1
%             -1, -1, 1
%             -1, -1, -1];        %Center of mass of each division
        

% F = [0, 0, 0];
% T = [0, 0, 0];
% for i=1:length(points(:,1))
%     r_p = r+points(i,:);
%     F = F + G*m1*m_TM/n/norm(r_p)^2.*r_p/norm(r_p);
%     T = T + cross(F,r_p);
% end
alpha = T/I;                    %Angular acceleration rad/s^2
a = F*m_TM;                     %Acceleration