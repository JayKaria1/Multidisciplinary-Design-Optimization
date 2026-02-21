% Constraints as per the homework file
function [c,ceq] = Constraints(X)

    Reference_Values;
    Initial_Geometry;
    
    global W_AminusW;
    global couplings;
    
    W_fuel = couplings.W_fuel;
    W_wing = couplings.W_wing;

    global initial;
    
    UCST = Ref.UCST .* (1 + X(1:6) * 0.07);
    LCST = Ref.LCST .* (1 + X(7:12) * 0.07);
    
    % Span (X13): lb = 0.86, ub = 1.28 → midpoint = 1.07, range/2 = 0.21
    Span = Ref.Span * (1.07 + X(13) * 0.21);
    
    % Area (X14): lb = 0.7, ub = 1.3 → midpoint = 1.0, range/2 = 0.3
    Area = Ref.Area * (1.0 + X(14) * 0.3);
    
    % TPR (X15): lb = 0.6, ub = 1.4 → midpoint = 1.0, range/2 = 0.4
    TPR = Ref.TR * (1.0 + X(15) * 0.4);
    
    % QC Sweep (X16): lb = 0.67, ub = 1.20 → midpoint = 0.935, range/2 = 0.265
    QC_Sweep = Ref.Sweep_25 * (0.935 + X(16) * 0.265);
    
    % M_cruise (X17): lb = 0.90, ub = 1.05 → midpoint = 0.975, range/2 = 0.075
    M_cruise = Ref.M_cruise * (0.975 + X(17) * 0.075);
    
    % H_cruise (X18): lb = 0.90, ub = 1.10 → midpoint = 1.0, range/2 = 0.10
    H_cruise = Ref.Alt * (1.0 + X(18) * 0.10);
    
    % spar_front (X19): lb = 0.83, ub = 1.11 → midpoint = 0.97, range/2 = 0.14
    spar_front = Ref.Spar_front * (0.97 + X(19) * 0.14);
    
    % spar_back (X20): lb = 0.95, ub = 1.03 → midpoint = 0.99, range/2 = 0.04
    spar_back = Ref.Spar_back * (0.99 + X(20) * 0.04);
        
    Current_MTOW= W_fuel+W_wing+W_AminusW;
    Init_MTOW = initial.W_fuel + initial.W_wing + W_AminusW;
    
    Reference_Wing_Loading = Init_MTOW/Ref.Area;
    Current_Wing_Loading = Current_MTOW/Area;
    
    c(1) = (Current_Wing_Loading - Reference_Wing_Loading)/Reference_Wing_Loading;   %Normalized the constraints w.r.t original Wing Loading
    
    Current_CO2 = W_fuel*3.16;
    Init_CO2 = initial.W_fuel * 3.16;
    
    c(2) = (Current_CO2 - Init_CO2)/Init_CO2; %Normialized w.r.t initial CO2 emission
    
    
    Current_Volume_Fuel = W_fuel/Ref.rho_fuel;
    Current_Fuel_Tank_Volume = FuelTank(TPR,Span,Area,QC_Sweep,RTA,TTA,Dihedral,spar_front,spar_back,UCST,LCST);
    
    c(3) = (Current_Volume_Fuel-Current_Fuel_Tank_Volume)/Ref.Volume_Fuel;
    
    ceq = [];

    fid = fopen('tolcheck.txt',"a");
    fprintf(fid,'%f\t%f\t%f\n',c(1),c(2),c(3));
    fclose(fid);

end




