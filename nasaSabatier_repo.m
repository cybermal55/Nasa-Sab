%define constants
A = 1.7e7; %L^4 mol^-4 s^-1, reaction rate prefactor
Ea = 80000; %m^3*Pa/mol, Activation energy CHECK UNITS!
Tnorm = 400 + 273.15; %degrees Kelvin, normal reactor temperature
T_reactor = 400 + 273.15; %K, temp inside reactor
P_reactor = 5; %bar, pressure inside reactor
R = 8.3145; %m^3*Pa/K*mol, gas constant
V_station = 1000000; %L, total pressurized volume of station
numAstronauts = 6; %six astronauts at station at any given time
CO2_perperson = 1.16*10^-5; %kg per second (500 L per day)
CO2_in = numAstronauts*CO2_perperson; %~6.9*10^-5
CO2_molarmass = .044009; %kg/mol
f_CO2 = 0.75; %conversion of CO2
%Stream1
xCO2_1 = 0.5;
xN2_1 = 0.5; %both at 25C and 2 bar
T_1 = 25 + 273.15; %K, temp of stream 1
P_1 = 2; %bar, pressure of stream 1
n_1 = 2*CO2_in/CO2_molarmass; %mol/s, kg/s of CO2 from astronauts/molar mass to be in mol/s, times 2 because 50% of the stream is N2
%Stream2
xH2_2 = 1.0; %0C and 5 bar
T_2 = 0 + 273.15; %K, temp of stream 2
P_2 = 5; %bar, pressure of stream 2
n_2 = 2*n_1; %mol/s, H2 stream is in stoichiometric quantity (4x as much as CO2)
%SR MUST CONVERT 75 % of the CO2 it receives to avoid harming astronauts

%%%ideal gas law assumption
%n_stream1 = CO@ of astronauts; %finding mols of stream 1
%n_stream2 = (P_2*V_station)/(R*T_2); %finding mols of stream 2 %4* amount of CO2
%%n_reactor = (P_reactor*V_R)/(R*T_reactor); %%%trL*bar/(K*molying to find V_R so cant do this

%pt. A
Fo = n_1 + n_2; %%inlet molar flow rate, SHOULD IT JUST BE n_1 BECAUSE THAT"S WHAT HAS THE CO2?
k = A*exp(-Ea/(R*Tnorm)); %L^4 mol^-4 s^-1, calculating k, using T_norm because this is at what temp reaction is taking place
vo = Fo*0.083145*T_reactor/P_reactor; %total V coming in reactor, L/s
%**WHICH V TO USE*
%V_1 = n_1*0.083145*T_1/P_1; %L units
%V_2 = n_2*0.083145*T_2/P_2; %L units
V_1 = n_1*0.083145*T_reactor/P_reactor; %L/s units
V_2 = n_2*0.083145*T_reactor/P_reactor; %L/s units
C_H2_0 = n_2*xH2_2/vo; 
C_CO2_0 = n_1*xCO2_1/vo; %finding concentration of CO2 in stream 1

C_CO2 = C_CO2_0*(1-f_CO2); 
C_H2 = C_H2_0 - 4*(C_CO2_0 - C_CO2); %rearranging C_H2 in terms of C_CO2 for rate law
negRate = k*C_CO2*(C_H2)^4; %% rate law !assuming stoich coeff -, can remove
V_R = (vo*C_CO2_0*f_CO2); 
V_R = V_R/negRate;

fprintf('Minimum volume of SR reactor is %f L\n', V_R)

%can now call on SR_odes that I solved for

C0 = [C_CO2_0; C_H2_0; 0; 0]; % initial concentrations in mol/L 
tspan = [0 3600]; % time span in seconds

[t, C] = ode45(@SR_odes, tspan, C0); %calling the function
%COMMENTING OUT 65-274 is giving me the right graphs
plot(t, C(:, 1)); %plotting the graph
xlabel('Time in seconds');
ylabel('Concentration of CO2 in mol/L');
title('Concentration of CO2 vs time');

%pt B continued
non_lethal_level = 0.03; %moles per liter
%uncomment if loop is preferred:
%V_R_guess = V_R; %iterate for part B
%V_R_decrease = 10000; %L I DID A ROOT FINDING METHOD BEFORE AND FOUND THE MINIMUM (commented below!) and it kept returning a value of 0L
%iteration_max = floor(V_R/ V_R_decrease); %(how many times it will iterate before volume is literally 0)
%for iterations = 1:iteration_max %number of iterations ... avoid max infinite loop
%
%[t, C] = ode45(@(t,C) SR_odes(t, C, n_2, V_R_guess), tspan, C0); %calling the function but changing ode based on the volume guess that we are iterating through
%
%C_CO2_out = C(end, 1) %final concentration of CO2
%disp(C_CO2_out)
% 
% if C_CO2_out <= non_lethal_level
%     fprintf('The minimum volume for the SR reactor to keep the CO2 levels as 0.03mol/L or below is %f L\n', V_R_guess); %print answer
%     break;
% else
%     V_R_guess = V_R_guess + V_R_decrease
%     C_CO2_out = C(end, 1); %final concentration of CO2
% end
% end
[~, idx] = min(abs(C_CO2-non_lethal_level));
Ctotal_out = C(end,:);
C_CO2_out = Ctotal_out(1); %final concentration from matrix
C_H2_out = Ctotal_out(2)

negRate_out = k*C_CO2_out*(C_H2_out)^4; %rate AT outlet at end of hour

V_R_belowtarget = (vo*C_CO2_out - non_lethal_level*vo)/negRate_out; %finding volume of SR to keep outlet conc at 0.03. check if reaction rate out is right!

fprintf('The minimum volume for the SR reactor to keep the CO2 levels at 0.03mol/L or below is %f L\n', V_R_belowtarget) %print answer

n_2_0 = 4*n_1;
% % 
figure; 
hold on;
multiples = [0.5, 1, 2, 4, 8]; %each factor of change that stoich. feed is multiplied by, that way loop can be used instead of doing what I did above which was tedious and messy
tspan = [0 3600]; %assuming t= an hour
C0 = [C_CO2_0; C_H2_0; 0; 0]; %initial conditions
colors = ['r', 'g', 'b', 'c', 'm']; %making it possible to plot everything in for loop

for i = 1:length(multiples)
    n_2 = multiples(i) * n_2_0;
    V_2_new = n_2 * R * Tnorm/P_reactor;
    C_H2_0_new = n_2 * xH2_2 / V_2_new;
    C0_new = [C_CO2_0; C_H2_0_new; 0; 0];

    [t,C] = ode45(@(t,C) SR_odes(t, C, n_2, V_R), tspan, C0_new);

    color = colors(i);
    plot(t,C(:,1),'color', color);
end

xlabel('Time in seconds');
ylabel('Concentration of CO2 in mol/L');
title('Concentration of CO2 over time with changing H2 inlet');
legend('0.5x', '1x', '2x', '4x', '8x');
grid on;

%Changing reactor temperature
figure; 
hold on;
temperature = [573, 623, 673, 723, 773]; %each diff T_reactor
tspan = [0 3600]; %assuming t= an hour
for i = 1:length(temperature)
    T_reactor = temperature(i);
    k = A*exp(-Ea/(R*T_reactor)); %L^4 mol^-4 s^-1, calculating k, using T_norm because this is at what temp reaction is taking place
    vo = Fo*0.083145*T_reactor/P_reactor;
    C_CO2_0 = n_1*xCO2_1/vo; %redefining initial concentration at new temps
    C_CO2 = C_CO2_0*(1-f_CO2);  %redefining concentrations at new temps
    C_H2_0 = n_2*xH2_2/vo; %same as above for C_H2
    C_H2 = C_H2_0 - 4*(C_CO2_0 - C_CO2);
    negRate = k*C_CO2*(C_H2)^4; %% rate law !assuming stoich coeff -, can remove
    V_R = (vo*C_CO2_0*f_CO2); 
    V_R = V_R/negRate;
    C0 = [C_CO2_0; C_H2_0; 0; 0];

    [t,C] = ode45(@(t,C) SR_odes(t, C, n_2, V_R), tspan, C0);

    color = colors(i);
    plot(t,C(:,1),'color', color);
end

xlabel('Time in seconds');
ylabel('Concentration of CO2 in mol/L');
title('Concentration of CO2 over time with changing reactor temperature in C');
legend('300', '350', '400', '450', '500');
grid on; %Note to the grader: I literally cannot for the life of me figure out why my concentration is staying constant despite all the changes, I have spent about 4.5 hours in total just trying to work out this issue alone and idk where I'm going wrong, I believe that this code would be right if the concentration was changing but for some reason it is stuck at 0.0605 (I think an issue with the ODE but I went line by line trying to remove the semi-colons and trouble shoot / adding in terms to see what would give me different results and everything either kept my values the same or errored:(


numAstronauts = numAstronauts + 2; %2 additional astronauts
%updating values using new number of astronauts
C_CO2_in = numAstronauts * CO2_perperson;
n_1 = 2*CO2_in/CO2_molarmass; %mol/s, kg/s of CO2 from astronauts/molar mass to be in mol/s, times 2 because 50% of the stream is N2
n_2 = 4*n_1;
Fo = n_1+n_2; %mol/s, this is technically FAo but for some reason I put it as vo
K = 0.05; %s^-1
V_MAR = 1; %L
V_1 = n_1*0.083145*T_reactor/P_reactor; %L units
V_2 = n_2*0.083145*T_reactor/P_reactor; %L units
C_H2_0 = n_2*xH2_2/V_2; %CHECK IF IDGL FOR V_1 and V_2 was correct!
C_CO2_0 = n_1*xCO2_1/V_1; %finding concentration of CO2 in stream 1
C_CO2 = C_CO2_0*(1-f_CO2); %outlet of SR for CO2
C_H2 = C_H2_0 - 4*(C_CO2_0 - C_CO2); %outlet of SR for H2
C_CH4 = C_CO2_0 - C_CO2; %assuming stoichiometric relationship for prod
C_H2O = 2*(C_CO2_0 - C_CO2); %assuming stoichiometric relationship for products

MAR1_inlet = Fo* (C_CO2 + C_H2 + C_CH4 + C_H2O); %total mol/s, Fo
F_CO2_0 = C_CO2 * Fo; %molar flow rate of JUST A into the MAR reactor
F_CO2_final = F_CO2_0* (1-f_CO2); %goal CO2 at the end of the MARs reactors

rateMAR = -K * C_CO2;
MAR_1_in_CO2 = C_CO2; %outlet of SR for CO2 is inlet of first MAR
count = 0;

while true
    count = count + 1;
    MAR_1_out_CO2 = MAR_1_in_CO2 * exp(-K*V_MAR/Fo); %this is the expression for C_CO2 in MARS because rate law is given as first order wrt CO2
    F_CO2 = Fo * MAR_1_out_CO2;
    if F_CO2 <= F_CO2_final
        fprintf('The number of MARs needed to react 75 percent of eight astronauts is %f \n', count);
        break;
    else
        MAR_1_in_CO2 = MAR_1_out_CO2;
    end
end



N = 100;
dz = V_MAR / N; %spatial step
u = vo / V_MAR; %superficial velocity (volumetric flow rate / volume)

CO2_mol_flow = CO2_in / CO2_molarmass; %mol/s
vo = Fo * R * T_reactor / P_reactor;   %L/s (ideal gas law)
C_CO2_in_molL = CO2_mol_flow / vo;     %mol/L

C_CO2_o = ones(N,1) * C_CO2_in_molL; %make sure this is in mol/L!

[t, C] = ode15s(@(t, C) PFR_MAR_ode(t, C, K, u, dz, N), [0 3600], C_CO2_o);

plot(t, C(:,end), 'r'); %plot outlet CO2 over time
xlabel('Time (s)');
ylabel('CO2 Concentration (mol/L)');
title('CO2 at MAR Reactor Outlet over Time');
%
F_CO2_spike = 0.75*n_1; %molar flow rate with spike
C_CO2_spike = F_CO2_spike / vo; %conc due to spike
C_CO2_MARs_in = C(end); %inlet of MARs is outlet of SR
N=100; %discretization points again
h = 0.1; %time step difference for Runge-Kutta
tspan = 0:h:3600; %assuming an hour for spike again
CO2_initial = C_CO2_spike * ones(N,1); %setting initial conditions for the CO2 spike for MARs reactors
CO2_initial(1) = C_CO2_MARs_in; %very first point is inlet of MARs which is outlet of SR

%manually implementing the 4th-order Runge-Kutta method
C_CO2_RK = zeros(length(tspan), N); %setting up a blank matrix
C_CO2_RK(1,:) = CO2_initial'; %setting all of the vals in col 1 every row (aka all initial CO2) equal to inlet MARs value

for idx = 1:length(tspan) - 1
    C = C_CO2_RK(idx,:)'; %setting up conc matrix using MARs concentration... transposing to pull using RK later
    K1 = -1*K*C; %defining change of k1 as of RK model
    K2 = -1*K*(C+0.5*h*K1); %k2 = f(x +h/2, y + h/2k1)
    K3 = -1*K*(C+0.5*h*K2); %k3 defined in problem sheet like above
    K4 = -1*K*(C+h*K3); %k4 defined in problem sheet like above

    C_CO2_RK(idx+1,:) = C_CO2_RK(idx,:) + (1/6)*h*(K1' + 2*K2' + 2*K3' +K4'); %formula for 4th order RK: y(t+h) = y(t) + (1/6) * (k1 + 2k2 + 2k3 + k4

    C_CO2_MARs_out = C_CO2_RK(end,N); %extracting final concentration based on last values of conc in the last (N) MARs reactor
end
    if C_CO2_MARs_out <= non_lethal_level
        fprintf('4th order RK method works. The number of reactor(s) needed to keep CO2 levels below 0.03 mol/L is %f \n', N);
    else 
        N = N+1;
    end
    % Plot CO2 Concentration Over Time for the last MAR
plot(tspan, C_CO2_RK, 'r'); % RK solution (Part D)
xlabel('Time in seconds');
ylabel('Concentration at exit of MARs');
title('Transient Behavior of CO2 Spikes in MAR (4th Order RK Method)');
grid on;
check = C(:,1);
check = check'
check2 = C_CO2_RK(:,end);
%Comparing results in parts B and D on the same plot

%C_interp = interp1(t, C(:,1), tspan, 'linear', 'extrap'); %extrapolate to make C(:,1) match the length of C_CO2_RK
C_partB = [0.019; C(:,1)];
last_val = C_partB(end);
C_extended = [C_partB(:,1); repmat(last_val, 36001 - length(C_partB), 1)];

figure;
newtime = reshape(tspan(1:100), [100, 1]);
plot(tspan, C_extended, 'r'); %RK solution (Part D)
hold on;
plot(tspan, C_CO2_RK(:,end), 'b'); %PDE solution (Part B)
xlabel('Time in seconds');
ylabel('Concentration at exit of MARs');
title('Transient Behavior of CO2 Spikes: PDE vs RK4 model');
legend('ODE45 solved (Part B)', 'RK4 solved (Part D)');
grid on;
hold off;

%tracking ALL species at outlet of last MAR

F_CO2_spike = 0.75 * n_1; % molar flow rate with spike
C_CO2_spike = F_CO2_spike / vo; % initial CO2 conc due to spike
C_CO2_MARs_in = C_CO2; % inlet of MARs is outlet of SR

N = 100; %discretization points again (spatial)
h = 0.1; %time step difference for Runge-Kutta
tspan = 0:h:3600; % assuming an hour for spike again

%IC for all species (assuming inlet from SR)
CO2_initial = C_CO2_spike * ones(N, 1); %CO2 conc profile
H2_initial = C_H2 * ones(N, 1);         %H2 conc profile (constant)
CH4_initial = C_CH4 * ones(N, 1);       %CH4 conc profile (constant)
H2O_initial = C_H2O * ones(N, 1);       %H2O conc profile (constant)

%storing the conc of all speicies over time
C_CO2_RK = zeros(length(tspan), N);
C_H2_RK = zeros(length(tspan), N);
C_CH4_RK = zeros(length(tspan), N);
C_H2O_RK = zeros(length(tspan), N);

%initial conditions
C_CO2_RK(1, :) = CO2_initial';
C_H2_RK(1, :) = H2_initial';
C_CH4_RK(1, :) = CH4_initial';
C_H2O_RK(1, :) = H2O_initial';

%RK4 Time loop
for idx = 1:length(tspan) - 1
    C_CO2 = C_CO2_RK(idx, :)'; %Transpose for K steps
    K1 = -K * C_CO2;
    K2 = -K * (C_CO2 + 0.5 * h * K1);
    K3 = -K * (C_CO2 + 0.5 * h * K2);
    K4 = -K * (C_CO2 + h * K3);
    C_CO2_new = C_CO2 + (h / 6) * (K1 + 2 * K2 + 2 * K3 + K4);

    %stoichiometry for other species at the outlet of the last MAR:
    delta_CO2 = C_CO2_new - C_CO2; % Change in CO2 concentration

    C_H2_new = C_H2_RK(idx, :)' - 4 * delta_CO2;
    C_CH4_new = C_CH4_RK(idx, :)' - delta_CO2;
    C_H2O_new = C_H2O_RK(idx, :)' - 2 * delta_CO2;

    %storing results for next time step
    C_CO2_RK(idx + 1, :) = C_CO2_new';
    C_H2_RK(idx + 1, :) = C_H2_new';
    C_CH4_RK(idx + 1, :) = C_CH4_new';
    C_H2O_RK(idx + 1, :) = C_H2O_new';
end

%extracting outlet concentrations (last position in the MARs)
C_CO2_outlet = C_CO2_RK(:, end);
C_H2_outlet = C_H2_RK(:, end);
C_CH4_outlet = C_CH4_RK(:, end);
C_H2O_outlet = C_H2O_RK(:, end);

%Plot ALL Species at outlet of the last MAR
figure;
hold on;
plot(tspan, C_CO2_outlet, 'r');
plot(tspan, C_H2_outlet, 'b');
plot(tspan, C_CH4_outlet, 'g');
plot(tspan, C_H2O_outlet, 'm');
xlabel('Time in seconds');
ylabel('Concentration of all active species in mol/L');
title('Outlet Concentrations over Time at Last MAR');
legend('CO2', 'H2', 'CH4', 'H2O');
grid on;

%- concentration of H2 (switched signs)
% RK method vs the PDE should converge
% Temperature graph
% varying H2 rate