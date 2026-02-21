%% Aerodynamic solver setting

% Wing planform geometry 

Root_chord = 5;
Tip_chord = 2;
Kink_chord = 3;
Inboard_span = 13;
Halfspan = 17;
Dihedral_angle = 5;
Root_twist = 5;
Kink_Twist = 3;
Tip_twist = -2;

Geom= Wing_planform(0.22,50,300,35,5,-2,5);
%                x    y     z   chord(m)    twist angle (deg) 

AC.Wing.Geom = [   Geom.Rx     ,    Geom.Ry ,  Geom.Rz    , Geom.Root_chord, 5;  % Root Leading Edge
           Geom.Kx, Geom.Ky, Geom.Kz , Kink_chord, Geom.Kink_Twist;         % Kink LE
           Geom.Tx    , Geom.Ty    , Geom.Tz     , Geom.Tip_chord , Tip_twist];         % Tip LE

%AC.Wing.Geom = WPG;
% Wing incidence angle (degree)
AC.Wing.inc  = 0;   
            
            
% Airfoil coefficients input matrix
%                    | ->     upper curve coeff.                <-|   | ->       lower curve coeff.       <-| 
AC.Wing.Airfoils   = [0.2171    0.3450    0.2975    0.2685    0.2893  -0.1299   -0.2388   -0.1635   -0.0476    0.0797;
                      0.2171    0.3450    0.2975    0.2685    0.2893  -0.1299   -0.2388   -0.1635   -0.0476    0.0797];
                  
AC.Wing.eta = [0;1];  % Spanwise location of the airfoil sections

% Viscous vs inviscid
AC.Visc  = 1;              % 0 for inviscid and 1 for viscous analysis
AC.Aero.MaxIterIndex = 150;    %Maximum number of Iteration for the
                                %convergence of viscous calculation
                                
                                
% Flight Condition
AC.Aero.V     = 250;            % flight speed (m/s)
AC.Aero.rho   = 0.4;         % air density  (kg/m3)
AC.Aero.alt   = 9600;             % flight altitude (m)
AC.Aero.Re    = 1.14e7*Geom.MAC;        % reynolds number (bqased on mean aerodynamic chord)
AC.Aero.M     = 0.2;           % flight Mach number 
% AC.Aero.CL    = 0.4;          % lift coefficient - comment this line to run the code for given alpha%
AC.Aero.Alpha = 2;             % angle of attack -  comment this line to run the code for given cl 


%% 
tic

Res = Q3D_solver(AC);

cd = Res.CDwing
toc