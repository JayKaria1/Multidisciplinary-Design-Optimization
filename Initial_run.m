function [W_fuel,W_wing,counter] = Initial_run(X, W_fuel_c, W_wing_c)
%For input X is normalized. 
%inputs are the inital guesses
%Defining the design variables X

    Reference_Values
    Initial_Geometry
    
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
    
    error = 1e-4;        %Tolerance
    W_fuel = 100000;
    W_wing = 20000;
    
    counter = 0;           %Iteration number

    global W_AminusW;
    global CD0;        
    
    global initial
    
    while abs((W_fuel - W_fuel_c))/W_fuel >error || abs(W_wing - W_wing_c)/W_wing > error
        
        if counter >0
            W_wing_c = W_wing;
            W_fuel_c = W_fuel;
        end
    
        [LD,CL,CD,Alpha] = Aero(Span,Area,TPR,QC_Sweep,UCST,LCST,W_wing_c,W_fuel_c,1);
        data = Loads(Span,Area,TPR,QC_Sweep,UCST,LCST,W_wing_c,W_fuel_c,H_cruise,1);
        W_wing = Structures(W_fuel_c,Span,Area,TPR,QC_Sweep,UCST,LCST,W_wing_c,spar_front,spar_back,1);
        W_fuel = Performance(M_cruise,H_cruise,W_wing_c,W_fuel_c,LD,1);
    
        counter = counter+1;
        % disp(counter)

        W_AminusW = Ref.MTOW - W_wing - W_fuel;
        CD0 = (CL/LD)-CD;
        
        initial.W_fuel = W_fuel;
        initial.W_wing = W_wing;

        
    end
    


end






