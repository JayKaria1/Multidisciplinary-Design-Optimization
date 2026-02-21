% Function File for performance. Created on 27th March 2025 by Arnav. All
% these formulas are found on pg 5. %Formula for W_fuel needed
% 0.938/W_fuelburn_ratio and brackets missing for prop_eff - 28/3 - Jay.
function W_fuel = Performance(M_cruise, alt_cruise, W_wing,W_fuel, LD,initial)
    Reference_Values;

    global W_AminusW;
    
    if initial==1
        MTOW = Ref.MTOW;
    else 
        MTOW = W_wing+W_fuel+W_AminusW;
    end
    [~,a,~,~] = atmosisa(alt_cruise);                                                      % ISA conditions at cruise altitude
    V = a*M_cruise;                                                                        % Cruise velocity in m/s
    prop_eff = exp((-((V - Ref.LRCruise)^2)/(2*(70^2)))-((alt_cruise-Ref.Alt)^2)/(2*(2500^2)));   % Propulsive efficiency formula as per MDO Documentation
    SFC = Ref.CT_bar/prop_eff;                                                                % Specific Fuel Consumption in N/Ns
    % Breguet Range Equation
    W_fuelburn_ratio = exp(Ref.DesignRange * SFC /(LD * V));                                        % Fuel burn in kg, where design range is constant. As per MDO documentation.
    W_fuel = (1 - (0.938/W_fuelburn_ratio))*MTOW;                                                % Fuel weight in kg, as per MDO Documentation. 
    CO2_emission = 3.16*W_fuel;                                                                     % CO2 emission in kg, as per MDO Documentation.