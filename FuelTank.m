function V_tank = FuelTank(TPR,Span,Area,QC_Sweep,RTA,TTA,Dihedral,spar_front,spar_back,UCST,LCST)
    %Takes into account fuel tank ends at 85%
    
    Extra_fuel_tank_volume = 0; %m3. Explaination in txt file for changes
    Geom = Wing_planform(TPR,Span,Area,QC_Sweep,RTA,TTA,Dihedral);

    X_cords = (0:0.01:1)';
    [Xu,Xl,c] = D_airfoil2(UCST,LCST,X_cords);
    
    Xu = Xu(:,2);
    Xl = Xl(:,2);
    
    sections = 100;
    spans = linspace(0,Geom.Halfspan*0.85,sections);
    
    chords = Geom.Root_chord-Geom.Root_chord*(1- TPR)*spans/Geom.Halfspan; % Formula JAY
    
    db = (Geom.Halfspan*0.85)/(sections+1);

    widths = (spar_back - spar_front)*chords; %Moved location since this needs spar locations in percentage. Below it is converted to indices

    spar_front = round(spar_front*100);
    spar_back = round(spar_back*100);
    
    frontspars = chords*(Xu(spar_front)-Xl(spar_front)); % Jay formula , denormalized
    rearspars = chords*(Xu(spar_back)-Xl(spar_back));
    % widths = (spar_back - spar_front)*chords; 
    
    areas = [];
    
    volume = 0;
    for i = 1:sections
        area = (frontspars(i)+rearspars(i))*widths(i)/2;
        volume = volume + area*db;
        areas(end+1) = area;
    end
    ftank = 0.93;
    
    V_tank = volume*2*ftank + Extra_fuel_tank_volume;
end