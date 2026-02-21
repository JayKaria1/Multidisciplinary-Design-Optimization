format longG
Initial_Geometry;
Reference_Values;
X0 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
[UCST0,LCST0,Span0,Area0,Sweep0,TPR0,M_cruise0,H_cruise0,spar_front0,spar_back0] = Denormalize(X0);

X_opt = load("Optimized_X.mat");
X_opt = X_opt.X_opt;

global W_AminusW;

[UCST_opt,LCST_opt,Span_opt,Area_opt,Sweep_opt,TPR_opt,M_Cruise_opt,H_cruise_opt,spar_front_opt,spar_back_opt] = Denormalize(X_opt);

Planform_init = Wing_planform(TPR0,Span0,Area0,Sweep0,RTA,TTA,Dihedral);
Planform_opt = Wing_planform(TPR_opt,Span_opt,Area_opt,Sweep_opt,RTA,TTA,Dihedral);

geom1 = Planform_init.WPG;
geom2 = Planform_opt.WPG;

wing_planform_plot(geom1,geom2);


[x1,y1] = CSTtoXY(UCST0,LCST0);
[x2,y2] = CSTtoXY(UCST_opt,LCST_opt);

% figure(Visible="off")
% plot(x2,y2);
% hold on;
% plot(x1,y1);
% exportgraphics(gcf,'plots/airfoil_comparison.pdf','ContentType','vector');
% hold off;
% 
% 
% 
% % ylim([-0.3 0.3]);
% % xlim([0 1.1]);
% 
% legend('Optimised airfoil','Initial airfoil');
% 
% plot_iso(UCST_opt,LCST_opt,Planform_init,Planform_opt)

%Plotting the lift drag distributions
% [W_fuel0,W_wing0,counter] = MDA_coordinator(X0, Ref.DFL, Ref.W_wing);
% [W_fuel_opt,W_wing_opt,counter_opt] = MDA_coordinator(X_opt,Ref.DFL,Ref.W_wing);

%Storing the numbers so that i dont have to run MDA every time. rerun it if
%unsure of values.
W_fuel0 = 1.178403125992397e+05;
W_wing0 = 3.52246e+04;
W_fuel_opt = 1.007607887545445e+05;
W_wing_opt = 2.24249e+04;

MTOW0 = W_fuel0 + W_wing0 + W_AminusW;
MTOW_opt = W_fuel_opt+W_wing_opt+W_AminusW;

[c0,ceq0] = Constraints_post(X0,W_fuel0,W_wing0,W_fuel0,W_wing0);
[c1,ceq1] = Constraints_post(X_opt,W_fuel_opt,W_wing_opt,W_fuel0,W_wing0);

W_fuel0 = Performance(M_cruise0, H_cruise0, W_wing0,W_fuel0, 16,0);
W_fuel_opt = Performance(M_Cruise_opt, H_cruise_opt, W_wing_opt,W_fuel_opt, 14.90675,0);

[Res0_loads] = Loads(Span0,Area0,TPR0,Sweep0,UCST0,LCST0,W_wing0,W_fuel0,H_cruise0,0);
[Res_opt_loads] = Loads(Span_opt,Area_opt,TPR_opt,Sweep_opt,UCST_opt,LCST_opt,W_wing_opt,W_fuel_opt,H_cruise_opt,0);
figure('Visible','off')
plot(Res0_loads.Wing.Yst./Planform_init.Halfspan,Res0_loads.Wing.ccl);
hold on;
plot(Res_opt_loads.Wing.Yst./Planform_opt.Halfspan,Res_opt_loads.Wing.ccl);
legend('Intial','Final')
xlabel('Spanwise Position');
ylabel('Spanwise Wing Lift')
exportgraphics(gcf,'plots/ccl_criticalcond.pdf','ContentType','vector')
hold off;


% [~,~,~,~,Res_inv0] = Aero_inviscid(Span0,Area0,TPR0,Sweep0,UCST0,LCST0,W_wing0,W_fuel0,0);
% [~,~,~,~,Res_inv_opt] = Aero_inviscid(Span_opt,Area_opt,TPR_opt,Sweep_opt,UCST_opt,LCST_opt,W_wing_opt,W_fuel_opt,0);
% 
% 
% 
% [LD0,~,~,~,Res0_aero] = Aero(Span0,Area0,TPR0,Sweep0,UCST0,LCST0,W_wing0,W_fuel0,0);
% [LD_opt,~,~,~,Res_opt_aero] = Aero(Span_opt,Area_opt,TPR_opt,Sweep_opt,UCST_opt,LCST_opt,W_wing_opt,W_fuel_opt,0);
% % figure('Visible','off')
% plot(Res0_aero.Wing.Yst./Planform_init.Halfspan,Res0_aero.Wing.ccl)
% hold on
% plot(Res_opt_aero.Wing.Yst./Planform_opt.Halfspan,Res_opt_aero.Wing.ccl)
% legend('Intial','Final')
% xlabel('Spanwise Position');
% ylabel('Spanwise Wing Lift')
% exportgraphics(gcf,'plots/ccl_designpoint.pdf','ContentType','vector')
% hold off;
% 
% chordlengthsinitial = [];
% % Your two different Y arrays
% yval_init = Res0_aero.Section.Y;  % e.g., Res0_aero.Section.Y
% yval_opt = Res_opt_aero.Section.Y;   % e.g., Res_opt_aero.Section.Y
% 
% % ===== Chord length for initial planform =====
% kink_pos_init = 0.35 * Planform_init.Halfspan;
% chordlengthsinitial = zeros(size(yval_init));
% 
% idx1_init = yval_init <= kink_pos_init;
% chordlengthsinitial(idx1_init) = Planform_init.Root_chord - ...
%     (Planform_init.Root_chord - Planform_init.Kink_chord) .* ...
%     (yval_init(idx1_init) / kink_pos_init);
% 
% idx2_init = yval_init > kink_pos_init;
% chordlengthsinitial(idx2_init) = Planform_init.Kink_chord - ...
%     (Planform_init.Kink_chord - Planform_init.Tip_chord) .* ...
%     ((yval_init(idx2_init) - kink_pos_init) / (Planform_init.Halfspan - kink_pos_init));
% 
% % ===== Chord length for optimized planform =====
% kink_pos_opt = 0.35 * Planform_opt.Halfspan;
% chordlengthsfinal = zeros(size(yval_opt));
% 
% idx1_opt = yval_opt <= kink_pos_opt;
% chordlengthsfinal(idx1_opt) = Planform_opt.Root_chord - ...
%     (Planform_opt.Root_chord - Planform_opt.Kink_chord) .* ...
%     (yval_opt(idx1_opt) / kink_pos_opt);
% 
% idx2_opt = yval_opt > kink_pos_opt;
% chordlengthsfinal(idx2_opt) = Planform_opt.Kink_chord - ...
%     (Planform_opt.Kink_chord - Planform_opt.Tip_chord) .* ...
%     ((yval_opt(idx2_opt) - kink_pos_opt) / (Planform_opt.Halfspan - kink_pos_opt));
% 
% 
% % spanwise wing drag coefficient ccd with separate curves for induced drag and profile + wave drag @ design point
% figure(Visible="off")
% plot(Res0_aero.Wing.Yst./Planform_init.Halfspan,Res0_aero.Wing.cdi.*Res0_aero.Wing.chord)
% hold on
% plot(Res0_aero.Section.Y./Planform_init.Halfspan,Res0_aero.Section.Cd .* transpose(chordlengthsinitial))
% hold on
% plot(Res_opt_aero.Wing.Yst./Planform_opt.Halfspan,Res_opt_aero.Wing.cdi.*Res_opt_aero.Wing.chord)
% hold on
% plot(Res_opt_aero.Section.Y./Planform_opt.Halfspan,Res_opt_aero.Section.Cd .* transpose(chordlengthsfinal))
% legend('Cdi intial', 'Cd initial', 'Cdi final', 'Cd final')
% xlabel('Spanwise Position')
% ylabel('Spanwise Wing Drag')
% exportgraphics(gcf,'plots/cd_designpoint.pdf','ContentType','vector')
% hold off

% Load the data robustly
% data = readmatrix('tolcheck.txt', 'Delimiter', '\t');
% 
% col1 = data(:,1);
% col2 = data(:,2);
% col3 = data(:,3);
% iterations = 1:length(col1);
% 
% % Create figure
% figure(Visible="off"); hold on;
% 
% % Define grey region from -Inf to 0.001
% ymin = min([col1; col2; col3]) - 0.001;
% ymax = 0.001;
% 
% x = [iterations, fliplr(iterations)];
% y = [ymax * ones(1, length(iterations)), ymin * ones(1, length(iterations))];
% 
% h = fill(x, y, [0.7 0.7 0.7]);
% set(h, 'FaceAlpha', 0.3, 'EdgeColor', 'none','HandleVisibility','off');
% 
% % Plot the columns
% plot(iterations, col1, 'r', 'LineWidth', 1, 'DisplayName', 'Wing Loading Constraint');
% plot(iterations, col2, 'g', 'LineWidth', 1, 'DisplayName', 'Emmissions Constraint');
% plot(iterations, col3, 'b', 'LineWidth', 1, 'DisplayName', 'Fuel Tank Volume Constraint');
% legend('Location','best')
% % Labels and grid
% 
% xlim([0,iterations]);
% xlabel('Iteration');
% ylabel('Value');
% % title('Column Values with Shaded y ≤ 0.001 Region');
% legend;
% grid on;
% 
% exportgraphics(gcf,'plots/Constraint_Plot.pdf','ContentType','vector')
% 
% hold off;



% xy_init = [0,0;
%     Planform_init.Ty,Planform_init.Halfspan;
%     Planform_init.Ty+Planform_init.Tip_chord,Planform_init.Halfspan;
%     Planform_init.Ky+Planform_init.Kink_chord,0.35*Planform_init.Halfspan;
%     Planform_init.Root_chord,0;
%     0,0    ];
% 
% xy_opt = [0,0;
%     Planform_opt.Ty,Planform_opt.Halfspan;
%     Planform_opt.Ty+Planform_opt.Tip_chord,Planform_opt.Halfspan;
%     Planform_opt.Ky+Planform_opt.Kink_chord,0.35*Planform_opt.Halfspan;
%     Planform_opt.Root_chord,0;
%     0,0    ];
% 
% rootx = x1.*Planform_init.Root_chord;
% rooty = y1.*Planform_init.Root_chord;
% rootz = transpose(repelem(0,length(x1)));
% kinkx = x1.*Planform_init.Kink_chord + Planform_init.Ky;
% kinky = y1.*Planform_init.Tip_chord;
% kinkz = transpose(repelem(Planform_init.Halfspan*0.35,length(x1)));
% tipx = x1.*Planform_init.Tip_chord+Planform_init.Ty;
% tipy = y1.*Planform_init.Tip_chord;
% tipz = transpose(repelem(Planform_init.Halfspan,length(x1)));
% figure
% thex = [rootx,kinkx,tipx];
% they = [rootz,kinkz,tipz];
% thez = [rooty,kinky,tipy];
% mesh(thex(1:50,:),they(1:50,:),thez(1:50,:),'EdgeColor','none','FaceColor','blue','FaceAlpha',0.5)
% axis equal
% hold on
% mesh(thex(51:end,:),they(51:end,:),thez(51:end,:),'EdgeColor','none','FaceColor','blue','FaceAlpha',0.5)
% 
% xoffset = 15;
% rootx = x2.*Planform_opt.Root_chord+xoffset;
% rooty = y2.*Planform_opt.Root_chord;
% rootz = transpose(repelem(0,length(x2)));
% kinkx = x2.*Planform_opt.Kink_chord + Planform_opt.Ky+xoffset;
% kinky = y2.*Planform_opt.Tip_chord;
% kinkz = transpose(repelem(Planform_opt.Halfspan*0.35,length(x2)));
% tipx = x2.*Planform_opt.Tip_chord+Planform_opt.Ty+xoffset;
% tipy = y2.*Planform_opt.Tip_chord;
% tipz = transpose(repelem(Planform_opt.Halfspan,length(x2)));
% hold on
% thex = [rootx,kinkx,tipx];
% they = [rootz,kinkz,tipz];
% thez = [rooty,kinky,tipy];
% mesh(thex(1:50,:),they(1:50,:),thez(1:50,:),'EdgeColor','none','FaceColor','red','FaceAlpha',0.5)
% axis equal
% hold on
% mesh(thex(51:end,:),they(51:end,:),thez(51:end,:),'EdgeColor','none','FaceColor','red','FaceAlpha',0.5)


function wing_planform_plot(geom1, geom2)
% Create figure
    figure(Visible="off"); hold on; grid on;
    
    % ---- Plot Geometry 1 ----
    x_LE1 = geom1(:,1);
    y_LE1 = geom1(:,2);
    x_TE1 = x_LE1 + geom1(:,4);
    plot([x_LE1; flipud(x_TE1)], [y_LE1; flipud(y_LE1)], 'k-', 'LineWidth', 2);
    scatter(x_LE1, y_LE1, 50, 'ro', 'filled');
    scatter(x_TE1, y_LE1, 50, 'rx', 'LineWidth', 1.5);
    
    % ---- Plot Geometry 2 ----
    x_LE2 = geom2(:,1);
    y_LE2 = geom2(:,2);
    x_TE2 = x_LE2 + geom2(:,4);
    plot([x_LE2; flipud(x_TE2)], [y_LE2; flipud(y_LE2)], 'b-', 'LineWidth', 2);
    scatter(x_LE2, y_LE2, 50, 'go', 'filled');
    scatter(x_TE2, y_LE2, 50, 'gx', 'LineWidth', 1.5);
    
    % Labels and formatting
    xlabel('x (Chordwise direction)');
    ylabel('y (Spanwise direction)');
    % title('Comparison of 2D Wing Planforms');
    
    legend('Initial Wing', '', '', ...
           'Optimized Wing', '', '', ...
           'Location', 'northwest');
    
    % Proper scaling for aircraft wings
    axis equal;
    
    % Add margins and set axis limits
    all_x = [x_LE1; x_TE1; x_LE2; x_TE2];
    all_y = [y_LE1; y_LE2];
    % x_margin = 0.1 * range(all_x);
    % y_margin = 0.05 * range(all_y);
    % 
    % xlim([min(all_x)-x_margin, max(all_x)+x_margin]);
    % ylim([min(all_y)-y_margin, max(all_y)+y_margin]);
    exportgraphics(gcf,'plots/2D_wing.pdf','ContentType','vector')
    hold off


end


function [UCST,LCST,Span,Area,QC_Sweep,TPR,M_cruise,H_cruise,spar_front,spar_back] = Denormalize(X)

    Reference_Values
    
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

end


function plot_iso(UCST_opt,LCST_opt,Planform_init,Planform_opt)

    [xAirfoil,yAirfoil] = CSTtoXY(UCST_opt,LCST_opt);
    
    % airfoilData = load('dsma523b.dat'); % [x, y] coordinates of the airfoil
    % xAirfoil = airfoilData(:, 1); % Airfoil x-coordinates (normalized)
    % yAirfoil = airfoilData(:, 2); % Airfoil y-coordinates (normalized)
    
    % W1 = load('WPG_final.mat');
    WPG1 = Planform_opt.WPG;
    
    % W2 = load("WPG_init.mat");
    WPG2 = Planform_init.WPG;
    
    X1 = []; Y1 = []; Z1 = [];
    X2 = []; Y2 = []; Z2 = [];
    
    % Plot the first WPG data
    figure(Visible="off");
    hold on;
    
    % Plot each airfoil from WPG1
    for i = 1:size(WPG1, 1)
        % Use each WPG point as the leading edge for the airfoil
        xLE = WPG1(i, 1);
        yLE = WPG1(i, 2);
        zLE = WPG1(i, 3);
        chord = WPG1(i, 4);
    
        % Scale airfoil by chord (no twist or translation)
        xScaled = xAirfoil * chord;
        yScaled = yAirfoil * chord;
    
        % Translate airfoil to the respective spanwise position
        X1 = [X1; xScaled + xLE];
        Y1 = [Y1; yScaled + yLE];
        Z1 = [Z1; zeros(size(xScaled)) + zLE];  % Flat at the z position of the leading edge
    
        % Plot each airfoil at root, kink, and tip for WPG1
        plot3(xScaled + xLE, yScaled + yLE, zeros(size(xScaled)) + zLE,'k-', 'LineWidth', 2);
    end
    
    % Highlight root, kink, and tip with lines connecting them for WPG1
    plot3(WPG1(:, 1), WPG1(:, 2), WPG1(:, 3), 'k-', 'LineWidth', 2, 'MarkerSize', 8);
    
    % Connect the leading edges between root, kink, and tip for WPG1 (solid lines)
    plot3([WPG1(1, 1), WPG1(2, 1)], [WPG1(1, 2), WPG1(2, 2)], [WPG1(1, 3), WPG1(2, 3)], 'k-', 'LineWidth', 2);
    plot3([WPG1(2, 1), WPG1(3, 1)], [WPG1(2, 2), WPG1(3, 2)], [WPG1(2, 3), WPG1(3, 3)], 'k-', 'LineWidth', 2);
    
    % Connect the trailing edges for WPG1 (solid lines)
    for i = 1:size(WPG1, 1)
        % Calculate trailing edge position (LE + chord length)
        xTE = WPG1(i, 1) + WPG1(i, 4);
        yTE = WPG1(i, 2);  % Same y as leading edge
        zTE = WPG1(i, 3);  % Same z as leading edge
        
        % Store trailing edge coordinates for WPG1
        WPG1_Trailing(i, :) = [xTE, yTE, zTE];
    end
    
    % Plot the trailing edge connections for WPG1 (solid lines)
    plot3([WPG1_Trailing(1, 1), WPG1_Trailing(2, 1)], [WPG1_Trailing(1, 2), WPG1_Trailing(2, 2)], ...
          [WPG1_Trailing(1, 3), WPG1_Trailing(2, 3)], 'k-', 'LineWidth', 2);
    plot3([WPG1_Trailing(2, 1), WPG1_Trailing(3, 1)], [WPG1_Trailing(2, 2), WPG1_Trailing(3, 2)], ...
          [WPG1_Trailing(2, 3), WPG1_Trailing(3, 3)], 'k-', 'LineWidth', 2);
    
    % Now plot the second WPG data (WPG2) over the first one
    
    % Plot each airfoil from WPG2
    for i = 1:size(WPG2, 1)
        % Use each WPG point as the leading edge for the airfoil
        xLE = WPG2(i, 1);
        yLE = WPG2(i, 2);
        zLE = WPG2(i, 3);
        chord = WPG2(i, 4);
    
        % Scale airfoil by chord (no twist or translation)
        xScaled = xAirfoil * chord;
        yScaled = yAirfoil * chord;
    
        % Translate airfoil to the respective spanwise position
        X2 = [X2; xScaled + xLE];
        Y2 = [Y2; yScaled + yLE];
        Z2 = [Z2; zeros(size(xScaled)) + zLE];  % Flat at the z position of the leading edge
    
        % Plot each airfoil at root, kink, and tip for WPG2
        plot3(xScaled + xLE, yScaled + yLE, zeros(size(xScaled)) + zLE, 'LineWidth', 2, 'Color', 'b');
    end
    
    % Highlight root, kink, and tip with lines connecting them for WPG2 (dotted lines)
    plot3(WPG2(:, 1), WPG2(:, 2), WPG2(:, 3), 'b-', 'LineWidth', 2, 'MarkerSize', 8);
    
    % Connect the leading edges between root, kink, and tip for WPG2 (dotted lines)
    plot3([WPG2(1, 1), WPG2(2, 1)], [WPG2(1, 2), WPG2(2, 2)], [WPG2(1, 3), WPG2(2, 3)], 'b-', 'LineWidth', 2);
    plot3([WPG2(2, 1), WPG2(3, 1)], [WPG2(2, 2), WPG2(3, 2)], [WPG2(2, 3), WPG2(3, 3)], 'b-', 'LineWidth', 2);
    
    % Connect the trailing edges for WPG2 (dotted lines)
    for i = 1:size(WPG2, 1)
        % Calculate trailing edge position (LE + chord length)
        xTE = WPG2(i, 1) + WPG2(i, 4);
        yTE = WPG2(i, 2);  % Same y as leading edge
        zTE = WPG2(i, 3);  % Same z as leading edge
        
        % Store trailing edge coordinates for WPG2
        WPG2_Trailing(i, :) = [xTE, yTE, zTE];
    end
    
    % Plot the trailing edge connections for WPG2 (dotted lines)
    plot3([WPG2_Trailing(1, 1), WPG2_Trailing(2, 1)], [WPG2_Trailing(1, 2), WPG2_Trailing(2, 2)], ...
          [WPG2_Trailing(1, 3), WPG2_Trailing(2, 3)], 'b-', 'LineWidth', 2);
    plot3([WPG2_Trailing(2, 1), WPG2_Trailing(3, 1)], [WPG2_Trailing(2, 2), WPG2_Trailing(3, 2)], ...
          [WPG2_Trailing(2, 3), WPG2_Trailing(3, 3)], 'b-', 'LineWidth', 2);
    
    % Set view and labels
        % Set view and axis properties
    view(45, 45);        % Isometric view
    axis equal;
    
    % Show only numeric ticks (no axis labels or grid)
    xlabel(''); ylabel(''); zlabel('');
    grid off;
    title('');

    % Enable ticks but remove tick labels
    xtickformat('%.0f');
    ytickformat('%.0f');
    ztickformat('%.0f');

    % Ensure ticks are visible
    set(gca, 'XColor', 'k', 'YColor', 'k', 'ZColor', 'k', ...
             'FontSize', 10, 'TickLength', [0.01 0.01]);

    % Export figure
    exportgraphics(gcf,'plots/wing_iso.pdf','ContentType','vector')
    hold off;


end


function [x,y] = CSTtoXY(UCST, LCST)
    % number of points for upper and lower side (so x2 in total)
    num_points = 500;

    % x-coordinates for which to calculate upper and lower y coords
    x_points = 0:1/(num_points-1):1;

    % D_airfoil2 expects vertical vector
    x_points = transpose(x_points); 

    % divide into separate upper and lower CST param arrays
    % CST_upper = CST(1:acf.CST_par_count/2);
    % CST_lower = CST(acf.CST_par_count/2+1:acf.CST_par_count);

    % run D_airfoil2
    [xyupper,xylower] = D_airfoil2(UCST,LCST,x_points);

    coords = [xyupper(length(xyupper):-1:1,:) % reverse
              xylower(2:length(xylower),:)]; % exclude first

    x = coords(:,1);
    y = coords(:,2);
end


% Constraints as per the homework file
function [c,ceq] = Constraints_post(X,W_fuel,W_wing,W_fuel0,W_wing0)

    Reference_Values;
    Initial_Geometry;
    
    global W_AminusW;
    % global couplings;
    % 
    % W_fuel = couplings.W_fuel;
    % W_wing = couplings.W_wing;

    % global initial;
    
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
    Init_MTOW = W_fuel0 + W_wing0 + W_AminusW;
    
    Reference_Wing_Loading = Init_MTOW/Ref.Area;
    Current_Wing_Loading = Current_MTOW/Area;
    
    c(1) = (Current_Wing_Loading - Reference_Wing_Loading)/Reference_Wing_Loading;   %Normalized the constraints w.r.t original Wing Loading
    
    Current_CO2 = W_fuel*3.16;
    Init_CO2 = W_fuel0 * 3.16;
    
    c(2) = (Current_CO2 - Init_CO2)/Init_CO2; %Normialized w.r.t initial CO2 emission
    
    
    Current_Volume_Fuel = W_fuel/Ref.rho_fuel;
    Current_Fuel_Tank_Volume = FuelTank(TPR,Span,Area,QC_Sweep,RTA,TTA,Dihedral,spar_front,spar_back,UCST,LCST);
    
    c(3) = (Current_Volume_Fuel-Current_Fuel_Tank_Volume)/Ref.Volume_Fuel;
    
    ceq = [];

    fid = fopen('tolcheck.txt',"a");
    fprintf(fid,'%f\t%f\t%f\n',c(1),c(2),c(3));
    fclose(fid);

end






