%Wing planform function file definition
function Geom = Wing_planform(TPR, Span, Area, QC_Sweep, Root_twist, Tip_twist, Dihedral_angle)

    Reference_Values;

    %Geom.Halfspan_Area = Area / 2;                           % Half of Wing Area
    Geom.Halfspan = Span / 2;                                % Half of Wing Span
    Inboard_span = 0.35 * 0.5*Ref.Span;                  % Spanwise Kink Distance
    %Inboard_Area = 0.35 * Geom.Halfspan_Area;                % Inboard Area
    %Outboard_Area = Geom.Halfspan_Area - Inboard_Area;       % Outboard Area
    TE_sweep_adjustment = 0.2 ;                          % Adjustment for Q3D, in degrees

    %Aspect Ratios for individual wing sections (Not needed, instead use wing area as a design variable)
    % Inboard_Aspect_Ratio = (Inboard_span^2) / Inboard_Area; % Inboard Aspect Ratio
    % Outboard_Aspect_Ratio = (Outboard_span^2) / Outboard_Area; % Outboard Aspect Ratio
   
    % Root and Tip Chords
    % Geom.Root_chord = (2*Area)/(Span*(1+TPR));  
    Geom.Root_chord = (Area+(Inboard_span^2 + (Inboard_span*(Geom.Halfspan-Inboard_span)))*tand(QC_Sweep-TE_sweep_adjustment))/(2*(Inboard_span+(Geom.Halfspan-Inboard_span)*(1+TPR)/2));
    Geom.MAC = 2*(Geom.Root_chord/3)*(1+TPR+TPR^2)/(1+TPR); % Mean Aerodynamic Chord
    Geom.Aspect_ratio = (Span^2)/Area;                      % Aspect Ratio

    %Leading Edge Sweep Calculation, returns in radians
    % Geom.LE_sweep = atan(tan(deg2rad(QC_Sweep))+0.25*(2*Geom.Root_chord/Span)*(1-TPR));
    Geom.LE_sweep = deg2rad(QC_Sweep);
    
   % Kink and Tip Chords
    Geom.Kink_chord = Geom.Root_chord - ( Inboard_span*(tan(Geom.LE_sweep)+tand(TE_sweep_adjustment)) ); 
    Geom.Tip_chord = Geom.Root_chord * TPR; 

    % Twist Angles
    Geom.Kink_Twist = Root_twist - (Root_twist - Tip_twist) * 0.35; 
    
    %Document: changed the way the output is because Q3D doesnt accept the
    %previous way for some reason 25/3 - Jay

    Geom.Rx = 0;
    Geom.Ry = 0;
    Geom.Rz = 0;

    Geom.Kx = Inboard_span * tan(Geom.LE_sweep);
    Geom.Ky =  Inboard_span * cosd(Dihedral_angle);
    Geom.Kz = Inboard_span * sind(Dihedral_angle);
    
    Geom.Tx =  Geom.Halfspan * tan(Geom.LE_sweep);
    Geom.Ty = Geom.Halfspan * cosd(Dihedral_angle);
    Geom.Tz =  Geom.Halfspan * sind(Dihedral_angle);

    % Wing Planform Geometry (WPG)

    Geom.WPG = [Geom.Rx, Geom.Ry, Geom.Rz, Geom.Root_chord, Root_twist;       % Root Leading Edge
                Geom.Kx, Geom.Ky, Geom.Kz, Geom.Kink_chord, Geom.Kink_Twist;  % Kink Leading Edge
                Geom.Tx, Geom.Ty, Geom.Tz, Geom.Tip_chord, Tip_twist];        % Tip Leading Edge

    % WPG = [            0                ,                    0               ,         0                          , Geom.Root_chord, Root_twist;         % Root Leading Edge
    %        Inboard_span * tan(Geom.LE_sweep), Inboard_span*(cosd(Dihedral_angle)), Inboard_span*(sind(Dihedral_angle)) , Geom.Kink_chord, Kink_Twist;         % Kink LE
    %        Geom.Halfspan * tan(Geom.LE_sweep)    , Geom.Halfspan*(cosd(Dihedral_angle))    , Geom.Halfspan*(sind(Dihedral_angle))     , Tip_chord , Tip_twist];         % Tip LE
end

