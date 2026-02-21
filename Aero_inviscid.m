%% Aerodynamic solver setting
% Used only for post-processing, inviscid analysis at design conditions
function [LD,CL,CD,Alpha,Res] = Aero_inviscid(Span,Area,TPR,QC_sweep,UCST,LCST,W_wing,W_fuel,initial)
    % Loading intitial values and reference values
    Reference_Values;
    Initial_Geometry;

    global W_AminusW;
    global CD0;                                   %Cd for aircraft fuselage

    if initial==1
        MTOW = Ref.MTOW;
    else 
        MTOW = W_wing+W_fuel+W_AminusW;
    end
    % Calling Wing planform function
    Geom = Wing_planform(TPR,Span,Area,QC_sweep,RTA,TTA,Dihedral);
    
    % AC.Wing naming as per Q3D Documentation (but it worked without it as well?)
    %Changed these as well. For some reason Geom.Kink_chord was giving an
    %error so I just replaced it with the formula. Messy but it works fine
    %now 25/3 -Jay
    AC.Wing.Geom = [            Geom.Rx   ,        Geom.Ry,    Geom.Rz     , Geom.Root_chord, 5;     % Root Leading Edge
                                Geom.Kx   ,        Geom.Ky,    Geom.Kz     , Geom.Root_chord - (0.15*Span * tan(Geom.LE_sweep)), Geom.Kink_Twist;         % Kink LE
                                Geom.Tx   ,        Geom.Ty,    Geom.Tz     , Geom.Tip_chord , -2];         % Tip LE;
    
    % Wing incidence angle (degree)
    AC.Wing.inc  = WIA;   
                
                
    % Airfoil coefficients input matrix
    %                    | ->     upper curve coeff.                <-|   | ->       lower curve coeff.       <-| 
    AC.Wing.Airfoils   = [Ref.UCST Ref.LCST;
                          Ref.UCST Ref.LCST;
                          Ref.UCST Ref.LCST];
    AC.Wing.eta = eta;             % Spanwise location of the airfoil sections at apex points which are root, kink and tip.
    
    % Viscous vs inviscid
    AC.Visc  = 0;                  % 0 for inviscid and 1 for viscous analysis
    AC.Aero.MaxIterIndex = 150;    % Maximum number of Iteration for the convergence of viscous calculation
          
    %Design Point (Cruise) ISA Conditions using the atmosisa package
    [T,a,P,rho] = atmosisa(Ref.Alt);           
    
    % Flight Condition definition
    L_des = sqrt(9.81*MTOW*(9.81*(MTOW-W_fuel)));                 % Formula in Pg 3 of MDO Document
    AC.Aero.V     = a*Ref.M_cruise;                               % flight speed (m/s)
    AC.Aero.rho   = rho;                                          % air density  (kg/m3)
    AC.Aero.alt   = Ref.Alt;                                      % flight altitude (m)
    AC.Aero.Re    = rho*AC.Aero.V*Geom.MAC/0.0000148;             % reynolds number (bqased on mean aerodynamic chord)
    AC.Aero.M     = AC.Aero.V/a;                                  % flight Mach number 
    AC.Aero.CL    = L_des/(0.5*rho*AC.Aero.V^2 * Area);           % lift coefficient - comment this line to run the code for given alpha%
    %AC.Aero.Alpha = 2;                                           % angle of attack -  comment this line to run the code for given cl 
    
    
    %% Timing the solver run
   
    
    Res = Q3D_solver(AC);         % Running the solver
    
    
    % Extracting the results (CD for induced drag case, hence CDi)
    CL = Res.CLwing;
    CD = Res.CDiwing;
    Alpha = Res.Alpha;
    
    if initial==1
        LD = 16;
    else
        LD = CL/(CD+CD0);
    end
end
