%% Aerodynamic solver setting
function [Res] = Loads(Span,Area,TPR,QC_sweep,UCST,LCST,W_wing,W_fuel,alt_cruise,initial)
    % Wing planform geometry 
    Reference_Values;
    Initial_Geometry;

    global W_AminusW;
    
    
    if initial==1
        MTOW = Ref.MTOW;
    else 
        MTOW = W_wing+W_fuel+W_AminusW;
    end

    Geom = Wing_planform(TPR,Span,Area,QC_sweep,RTA,TTA,Dihedral);
    
    % AC.Wing naming as per Q3D Documentation
    %Changed these as well. For some reason Geom.Kink_chord was giving an
    %error so I just replaced it with the formula. Messy but it works fine
    %now 25/3 -Jay
    AC.Wing.Geom = [            Geom.Rx   ,       Geom.Ry   ,     Geom.Rz   ,   Geom.Root_chord, 5;     % Root Leading Edge
                                Geom.Kx   ,       Geom.Ky   ,     Geom.Kz   ,   Geom.Root_chord - (0.15*Span * tan(Geom.LE_sweep)), Geom.Kink_Twist;         % Kink LE
                                Geom.Tx   ,       Geom.Ty   ,     Geom.Tz   ,   Geom.Tip_chord , -2];         % Tip LE;
    % Wing incidence angle (degree)
    AC.Wing.inc  = WIA;   
                
                
    % Airfoil coefficients input matrix
    %                    | ->     upper curve coeff.                <-|   | ->       lower curve coeff.       <-| 
    AC.Wing.Airfoils   = [Ref.UCST Ref.LCST;
                          Ref.UCST Ref.LCST;
                          Ref.UCST Ref.LCST];

    AC.Wing.eta = [0;0.35;1];  % Spanwise location of the airfoil sections
    
    % Viscous vs inviscid
    AC.Visc  = 0;                  % 0 for inviscid and 1 for viscous analysis
    AC.Aero.MaxIterIndex = 150;    % Maximum number of Iteration for the
                                   % convergence of viscous calculation
                                    
    [~,a,~,rho] = atmosisa(Ref.Alt);                               
    % Flight Condition
    L_max = 9.81*MTOW*2.5;
    AC.Aero.V     = Ref.MMO*a;                              % flight speed (m/s), where load factor is 2.5
    AC.Aero.rho   = rho;                                     % air density  (kg/m3)
    AC.Aero.alt   = alt_cruise;                              % flight altitude (m)
    AC.Aero.Re    = rho*AC.Aero.V*Geom.MAC/0.0000148;        % reynolds number (based on mean aerodynamic chord)
    AC.Aero.M     = Ref.MMO;                                     % AC.Aero.V/a;           % Maximum operating Mach number 
    AC.Aero.CL    = L_max/(0.5*rho*AC.Aero.V^2 * Area);      % lift coefficient - comment this line to run the code for given alpha%
    %AC.Aero.Alpha = 2;             % angle of attack -  comment this line to run the code for given cl 
    
    
    %% 
    
    
    Res = Q3D_solver(AC);
   
    
    ccl = Res.Wing.ccl(:);
    y = Res.Wing.Yst(:);
    cm4 = Res.Wing.cm_c4(:);
    chords = Res.Wing.chord(:);
    y_norm = y/(0.5*Span); %
    yi = linspace(0,1,16);

    cclo = interp1(y_norm,ccl,yi,"spline",'extrap');
    cm4o =interp1(y_norm,cm4,yi,"spline",'extrap'); 
    chordso = interp1(y_norm,chords,yi,"spline",'extrap');

    q = 0.5*AC.Aero.rho*AC.Aero.V^2;    

    L = [];
    M = [];

    for i = 1:length(yi)
        
        l = cclo(i)*q;
        m = chordso(i)*Geom.MAC*cm4o(i)*q;
        L(end+1) = l;
        M(end+1) = m;
    end
    data = [yi', L', M'];

    output_file = 'DC10.load';

    % Open the file for writing
    fid = fopen(output_file, 'w');
    
    % Write the header
    %fprintf(fid, 'spanwise position Lift Pitching moment\n');
    
    % Write the data in the required format
    fprintf(fid, '%.4f %.4e %.4e\n', data');
    
    % Close the file
    fclose(fid);

end
