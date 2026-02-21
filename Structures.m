function W_wing = Structures(W_fuel,Span,Area,TPR,QC_sweep,UCST,LCST,W_wing,Spar_front, Spar_rear,initial)
    Initial_Geometry;
    Reference_Values;
    Airfoil_to_EMWET(UCST,LCST);

    global W_AminusW;
    if initial==1
        MTOW = Ref.MTOW;
    else
        MTOW = W_wing+W_fuel+W_AminusW;
    end
    
    Geom = Wing_planform(TPR,Span,Area,QC_sweep,RTA,TTA,Dihedral);
    
    namefile    =    char('DC10');
    MTOW        =    MTOW;         %[kg]
    %MZF         =    W_AminusW + W_wing - Ref.DP;         %[kg] Need
    %W_AminusW which is a constant. commented for now since we need to do
    %initial run
    MZF = MTOW-W_fuel;
    % if initial==1
    %     MZF = Ref.ZFW;
    % end


    nz_max      =    2.5;   % Load Factor
    span        =    Span;            %[m]
    root_chord =    Geom.Root_chord;           %[m]
    taper       =    TPR;          
    sweep_le    =    rad2deg(Geom.LE_sweep);             %[deg]
    spar_front  =    Spar_front;
    spar_rear   =    Spar_rear;
    ftank_start =    0;                 % Starts at root chord as per pg 6 of the homework pdf
    ftank_end   =    0.85;              % Till 85 pc of the wing span
    eng_num     =    1;
    eng_ypos    =    0.3;
    eng_mass    =    Ref.Engine_mass;         %kg
    E_al        =    7E10;       %N/m2
    rho_al      =    2800;         %kg/m3
    Ft_al       =    295*10^6;        %N/m2 Tension Yield Stress
    Fc_al       =    295*10^6;        %N/m2 % Compression Yield Stress
    pitch_rib   =    0.5;          %[m]
    eff_factor  =    0.96;             % Depend on the stringer type
    Airfoil     =    'withcomb135_EMWET'; % update based on cst not airfoil
    section_num =    3;
    airfoil_num =    3;
    wing_surf   =    Area;
    
    fid = fopen( 'DC10.init','wt');
    fprintf(fid, '%g %g \n',MTOW,MZF);
    fprintf(fid, '%g \n',nz_max);
    
    fprintf(fid, '%g %g %g %g \n',wing_surf,span,section_num,airfoil_num);
    
    fprintf(fid, '0 %s \n',Airfoil);
    fprintf(fid, '0.35 %s \n',Airfoil);
    fprintf(fid, '1 %s \n',Airfoil);
    fprintf(fid, '%g %g %g %g %g %g \n',Geom.Root_chord,Geom.Rx,Geom.Ry,Geom.Rz,spar_front,spar_rear);
    fprintf(fid, '%g %g %g %g %g %g \n',Geom.Kink_chord,Geom.Kx,Geom.Ky,Geom.Kz,spar_front,spar_rear);
    fprintf(fid, '%g %g %g %g %g %g \n',Geom.Tip_chord,Geom.Tx,Geom.Ty,Geom.Tz,spar_front,spar_rear);
    
    fprintf(fid, '%g %g \n',ftank_start,ftank_end);
    
    fprintf(fid, '%g \n', eng_num);
    fprintf(fid, '%g  %g \n', eng_ypos,eng_mass);
    
    fprintf(fid, '%g %g %g %g \n',E_al,rho_al,Ft_al,Fc_al);
    fprintf(fid, '%g %g %g %g \n',E_al,rho_al,Ft_al,Fc_al);
    fprintf(fid, '%g %g %g %g \n',E_al,rho_al,Ft_al,Fc_al);
    fprintf(fid, '%g %g %g %g \n',E_al,rho_al,Ft_al,Fc_al);
    
    fprintf(fid,'%g %g \n',eff_factor,pitch_rib);
    fprintf(fid,'0 \n');
    fclose(fid);
    
    
    
    EMWET DC10;

    resultfile = fopen('DC10.weight','r');
    line = fgetl(resultfile);
    W_wing_array = textscan(line,"Wing total weight(kg)%f");
    W_wing= (W_wing_array{1});
    
end