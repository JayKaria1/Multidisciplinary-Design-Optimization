% This file contains all the values in the excel file provided to us by the professors. This file should not be tinkered with.

Ref.Thrust = 236*(10^3);      %N (Static Thrust)
Ref.Ramp_Weight = 265000;     %kg 
Ref.MTOW = 263636;            %kg Maximum Take-off Weight
Ref.MLW = 182798;             %kg Maximum Landing Weight
Ref.ZFW = 166922;             %kg Zero Fuel Weight
Ref.MPL = 48330;              %kg Maximum Payload
Ref.MFP = 32000;              %kg Maximum Payload at Maximum Fuel
Ref.DP = 26315;               %kg Design Payload
Ref.DFL = 115957;             %kg Design Fuel Load
Ref.OEW = 121364;             %kg Operating Empty Weight
Ref.rho_fuel = 0.81715*10^3;       %kg/litre Fuel Density given in MDO Document page number 6
Ref.ftank = 0.93;                          %Fuel Tank Volume Fraction given in MDO Document page number 6
Ref.FuelCapacity = 138165*Ref.rho_fuel;    %maximum mass Standard Version (kg). not used anywhere
Ref.Engine_mass = 4100;       %kg Engine Mass
Ref.W_wing = 0.12*Ref.MTOW ; %Initial guess for Wing Weight for Initial run based on engineering judgement


Ref.FusLength = 51.97;        % m Fuselage Length
Ref.FusHeight = 6.02;         % m Fuselage Height
Ref.FusWeidth = 6.02;         % m Fuselage Width 
Ref.Finess = 8.63;            % Fuselage Finess Ratio 

Ref.Area = 367.70;                    % This should be parameterized as 2 trapezoids and done somewhere else
Ref.Wing_loading = Ref.MTOW/Ref.Area; % Reference wing loading factor , goes into the constraints file
Ref.Span = 50.40;             % m Wingspan
Ref.MAC = 8.59;               % m Mean Aerodynamic Chord
Ref.AR = 6.91;                % Aspect Ratio
Ref.TR = 0.22;                % Taper Ratio
Ref.tc = 0.11;                % Average Thickness-to-Chord Ratio
Ref.Sweep_25 = 35;            % deg Quarter Chord Sweep
Ref.CD0 = 0; %Initial Guess for Fuselage Drag
Ref.LD = 16; %Initial L/D Ratio from assignment

Ref.V2 = 94.6577696;          % m/s V2 Speed (Rotation at TO)
Ref.Vapp = 75.6233268;        % m/s Approach Speed
Ref.MMO = 0.88;               % Maximum Operating Mach Number
Ref.CLmaxTO = 1.85;           % Maximum CL at TakeOff
Ref.CLmaxL = 2.35;            % Maximum CL at Landing 
Ref.M_cruise = 0.76;%Reverted to mach 0.78Try M = 0.76 since airfoil change may not be resolving in q3d when in fmincon  % 0.835 in ref excel, but Q3D can't handle that % Cruise mach number changed manually because Q3D doesnt work with given value
Ref.MaxSpeed = 272.655;       % m/s Maximum Cruise Speed 
Ref.LRCruise = 252.078;       % m/s Cruise Speed for Long Range
Ref.Alt = 31000 * 0.3048;     % m Cruise Altitude

%Units in m[m]. Was converted wrongly before 28/3 - Jay
Ref.RangeMPL = 4000*1852;          % in m Range and Max Payload
Ref.DesignRange  = 5373*1852;  % in m Design Range 
Ref.RangeMFP = 6400*1852;        % in m Range at Maximum Fuel + Payload

%Performance Module Specific Reference Values
Ref.CT_bar = 1.8639*(10^-4);           % CT bar value in N/NS, as those are the units used in Brueget Range Equation.
Ref.Spar_front = 0.1; %Just changed these to see if fuel tank constraint matches
Ref.Spar_back = 0.65;
Ref.CO2 = 3.16*Ref.DFL;
Ref.f_tank = 0.93;                     %From Page 6 in Homework

Ref.Volume_Fuel = Ref.DFL/Ref.rho_fuel;

%Airfoil Paramtrization              % UB           LB
Ref.LCST1 = -0.225404261885337;          % -0.2        -0.24
Ref.LCST2 = -0.163420361373538;          % -1.47       -1.8
Ref.LCST3 = -0.0470228855002990;         % -0.043      -0.051
Ref.LCST4 = -0.477077492806743;          % -0.43       -0.51
Ref.LCST5 = 0.0734577258998625;          % 0.08         0.066
Ref.LCST6 = 0.325538118792759;           % 0.36         0.28
%Upper Airfoil Coefficients
Ref.UCST1 = 0.233729303488509;           % 0.25         0.21
Ref.UCST2 = 0.0795985746292206;          % 0.087        0.071
Ref.UCST3 = 0.268260966221277;           % 0.29         0.23
Ref.UCST4 = 0.0887220843087906;          % 0.097        0.08
Ref.UCST5 = 0.278863675812943;           % 0.3          0.24
Ref.UCST6 = 0.381120210230920;           % 0.42         0.34             % These were done for 10% bounds. Not suitable due to instersection. Bounds must be 7%

Ref.LCST = [Ref.LCST1 Ref.LCST2 Ref.LCST3 Ref.LCST4 Ref.LCST5 Ref.LCST6];                %Lower Airfoil Coefficients
Ref.UCST = [Ref.UCST1 Ref.UCST2 Ref.UCST3 Ref.UCST4 Ref.UCST5 Ref.UCST6];                %Upper Airfoil Coefficients
