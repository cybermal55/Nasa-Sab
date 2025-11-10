function dCdt = SR_odes(t, C, n_1, n_2, V_R, T_reactor)
    % Constants
    A = 6.2e-8;            % L^4 mol^-4 s^-1
    Ea = 80000;            % J/mol (adjust units if needed!)
    R = 8.3145;            % J/mol/K
    vo = n_1 + n_2;        % total molar flow rate

    % Initial concentrations
    xCO2_1 = 0.75;
    xH2_2 = 1.0;

    C_CO2_0 = n_1 * xCO2_1 / V_R; % mol/L
    C_H2_0 = n_2 * xH2_2 / V_R;   % mol/L

    % Concentrations from state vector
    C_CO2 = C(1);
    C_H2  = C(2);
    C_CH4 = C(3);
    C_H2O = C(4);

    % Temperature-dependent rate constant
    k = A * exp(-Ea / (R * T_reactor)); % L^4/mol^4/s

    % Reaction rate
    r = k * C_CO2 * C_H2^4;

    % Material balances (assume ideal plug flow)
    dC_CO2 = -r + (C_CO2_0 - C_CO2) * (vo / V_R);
    dC_H2  = -4 * r + (C_H2_0  - C_H2) * (vo / V_R);
    dC_CH4 =  r + (-C_CH4) * (vo / V_R);
    dC_H2O = 2*r + (-C_H2O) * (vo / V_R);

    dCdt = [dC_CO2; dC_H2; dC_CH4; dC_H2O];
end



% function dCdt = SR_odes(t,C, n_2, V_R)
% %%part B
% %Defining constants
%     A = 6.2e-8; %L^4 mol^-4 s^-1, reaction rate prefactor
%     Ea = 80000; %m^3*Pa/mol, Activation energy CHECK UNITS!
%     Tnorm = 400 + 273.15; %degrees Kelvin, normal reactor temperature
%     T_reactor = 400 + 273.15; %K, temp inside reactor
%     P_reactor = 5; %bar, pressure inside reactor
%     R = 8.3145; %m^3*Pa/K*mol, gas constant
%     V_station = 1000000; %L, total pressurized volume of station
%     numAstronauts = 6; %six astronauts at station at any given time
%     CO2_perperson = 1.16*10^-5; %kg per second (500 L per day)
%     CO2_in = numAstronauts*CO2_perperson; %~6.9*10^-5
%     CO2_molarmass = .044009; %kg/mol
%     f_CO2 = 0.75; %conversion of CO2
% %Stream1
%     xCO2_1 = 0.75;
%     T_1 = 25 + 273.15; %K, temp of stream 1
%     P_1 = 2; %bar, pressure of stream 1
%     n_1 = 2*CO2_in/CO2_molarmass; %mol/s, kg/s of CO2 from astronauts/molar mass to be in mol/s, times 2 because 50% of the stream is N2
%     V_R = vo*0.083145*T_reactor/P_reactor; %total V coming in reactor
%     %Stream2
%     xH2_2 = 1.0; %0C and 5 bar
%     T_2 = 0 + 273.15; %K, temp of stream 2
%     P_2 = 5; %bar, pressure of stream 2
%     n_2 = 4*n_1; %mol/s, H2 stream is in stoichiometric quantity (4x as much as CO2)
% %SR MUST CONVERT 75 % of the CO2 it receives to avoid harming astronauts
%     vo = n_1 + n_2; %%inlet flow rate, SHOULD IT JUST BE n_1 BECAUSE THAT"S WHAT HAS THE CO2?
%     k = A*exp(-Ea/(R*Tnorm));
%     V_1 = n_1*0.083145*T_1/P_1; %CHANGE THIS BASED ON OUR V_1
% 
%     C_CO2_0 = n_1*xCO2_1/V_R; %finding concentration of CO2 in stream 1
%     C_H2_0 = n_2*xH2_2/V_R
% %%%%%%%%%%%%%%%%%%%%%%%%
% dCdt = zeros(4,1); %empty matrix
% 
% C_CO2 = C(1); %getting the concentration out of the matrix
% C_H2 = C(2);
% C_CH4 = C(3);
% C_H2O = C(4);
% 
% negRate = -1*k*C_CO2*(C_H2)^4; %% rate law
% 
% %defining the diffeqs
% % dC_CO2/dt = -k*C_CO2*(C_H2)^4-C_CO2_0*n_1 - C_CO2*vo %% should I use vo? n_1?
% % dC_H2/dt= -4*k*C_CO28(C_H2)^4-C_H2_0*n_2 - C_H2*vo
% % dC_CH4/dt = k*C_CO2*(C_H2)^4-C_CH4*vo
% % dC_H2O/dt = 2*k*C_CO2*(C_H2)^4-C_H2O*vo
% 
% dCdt(1) = negRate  + (C_CO2_0 - C_CO2)*(vo/V_R); %CHECK SIGNS!rate of CO2 going to products, - CO2 coming in, + CO2 leaving (should signs be swapped?)
% dCdt(2) = 4*negRate+ (C_H2_0  - C_H2)*(vo/V_R);  %rate of H2 going to products, - H2 coming in, + H2 leaving (should signs be swapped?)
% dCdt(3) = -1*negRate + C_CH4*(vo/V_R); %rate of CH4 being formed, - CH4 leaving (should signs be swapped?)
% dCdt(4) = -2*negRate + C_H2O*(vo/V_R); %rate of H2O being formed, - H2O leaving (should signs be swapped?)
% 
% dCdt = [dCdt(1); dCdt(2); dCdt(3); dCdt(4)];
% end
% 
