function dCdt = PFR_MAR_ode(t, C, K, u, dz, N)
    dCdt = zeros(N, 1); %preallocate derivative vector

    %boundary condition at inlet
    dCdt(1) = -u * (C(1) - C(2)) / dz - K * C(1);

    % Interior nodes (i = 2 to N-1)
    for i = 2:N-1
        dCdz = (C(i) - C(i+1)) / dz;  %backward difference
        dCdt(i) = -u * dCdz - K * C(i);
    end

    %outlet node(node N)
    dCdt(N) = -u * (C(N) - C(N-1)) / dz - K * C(N);
end


% function dCdt = PFR_MAR_ode(t, C, u, dz, K, MAR1_inlet)
%     %defining constants since this is a new function
%     V_MAR = 1; %L
%     K = 0.05; %s^-1
%     numAstronauts = 8; %number of atronauts
%     CO2_perperson = 1.16*10^-5; %kg per second (500 L per day)
%     C_CO2_in = numAstronauts * CO2_perperson; %CO2 coming into SR reactor
%     C_CO2_spike = 0.0194; %mol/L, F_CO2_spike / vo but I just ran it without ; in other script to get value so I wouldn't have to paste over 4 diff lines of code
%     non_lethal_level = 0.03; %mol/L
%     N = 100; %discretization of points
%     dz = V_MAR / N; %step size 
% 
%     tspan = [0 3600]; %assuming an hour
%     dt = 1; %time step
% 
%     MARs_quantity = 1; %initial number of MARs
%     C_CO2_outlet = C_CO2_spike;
% 
%     dCdt = zeros(N, 1); %initial empty matrix 
%     dCdt(1) = -1*K*C(1) - (C_CO2_in/dz)*(C(1)-C(2)); %%using guidance from https://www.youtube.com/watch?v=ZxfmmVAE2l4 video!!!
% 
%     for i = 2:N-1
%         dCdt(i) = -1*K*C(i) - (C_CO2_in/dz)*(C(i)-C(i+1)); %repeat iteration 
%     end
% % 
% %     dCdt(N) = -1*K*C(N)-(C_CO2_in/dz)*(C(N)-0);%end discretization
% % end
% % 
% %     N = length(C); %defining spacial descretization points
% %     dCdt = zeros(N, 1); %creating an empty matrix to reference back yo
% %     %C(1) = MAR1_inlet %initializing matrix by putting first term as inlet of MARs reactor
% %     dCdt(1) = u * (MAR1_inlet - C(1)) / deltaz - K * C(1); %upper bound/initial condition
% % 
% %     for i = 2:N
% %         dCdz = (C(i) - C(i-1)) / deltaz;
% %         % dCdz(i) = (C(i)-C(i-1))/deltaz; 
% %         dCdt(i) = -u * dCdz - K * C(i); %discretization of nodes
% % 
% %     end 
% %     % dCdt(1)= 0 %assuming initial conc at very beginning of MArs is close enough to 0
% % end
