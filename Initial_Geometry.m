%CST Coefficients for Whitcomb Airfoil

%Initial Values and constant values
Initial_span = 50.4;                   % Initial Wing Span Constant
RTA = 3;                               % Root Twist Angle Constant
TTA = -2;                              % Tip Twist Angle Constant
WIA = 2;                               % Wing Incidence Angle Constant
Initial_MAC = 8.59;                    % Mean Aerodynamic Chord 
Initial_TPR = 0.22;                    % Taper Ratio 
Initial_area = 367.70;                 % Wing Area
Inboard_span = 0.35*0.5*Initial_span;  % Inboard Spanwise Kink Distance Constant
Initial_QCsweep = 35;                  % Initial Quarter Chord Sweep
Dihedral = 5;                          % Dihedral Angle Constant


eta = [0;0.35;1]; %Position of airfoils 

% Caaling the Wing Planform Function
Geom = Wing_planform(Initial_TPR, Initial_span, Initial_area , Initial_QCsweep, RTA, TTA, Dihedral);

% Plotting the Wing Planform

data = Geom.WPG;
% Extract leading edge (LE) coordinates
x_LE = data(:,1);   % x-coordinates (first column)
y_LE = data(:,2);   % y-coordinates (second column)

% Extract chord lengths
c = data(:,4);      % chord length (fourth column)

% Compute trailing edge (TE) coordinates
x_TE = x_LE + c;  % TE is at LE minus chord length

%Plot the wing shape
% figure; hold on; grid on;
% plot([x_LE; flipud(x_TE)], [y_LE; flipud(y_LE)], 'k-', 'LineWidth', 2); % Wing outline
% scatter(x_LE, y_LE, 50, 'ro', 'filled'); % Mark LE points
% scatter(x_TE, y_LE, 50, 'bo', 'filled'); % Mark TE points
% 
% % Labels and formatting
% xlabel('x (Chordwise direction)');
% ylabel('y (Spanwise direction)');
% title('2D Wing Planform');
% legend('Wing Outline', 'Leading Edge', 'Trailing Edge', 'Location', 'Best');
% axis equal; % Maintain aspect ratio