function Airfoil_to_EMWET(UCST,LCST)

% Read the airfoil data from the input file
% input_file = 'withcomb135.dat';      % Original file
% output_file = 'withcomb135_EMWET.dat'; % Reformatted file

% Initial_Geometry;

% Load data (assuming two-column format: x-coordinates and y-coordinates)
    num_points = 50;

    % x-coordinates for which to calculate upper and lower y coords
    x_points = 0:1/(num_points-1):1;

    % D_airfoil2 expects vertical vector
    x_points = transpose(x_points); 

    % divide into separate upper and lower CST param arrays
    % CST_upper = Y.CST(1:acf.CST_par_count/2);
    % CST_lower = Y.CST(acf.CST_par_count/2+1:acf.CST_par_count);

    % run D_airfoil2
    [xyupper,xylower] = D_airfoil2(UCST,LCST,x_points);


    %
    % write result to .dat file
    %

    % assemble upper and lower coords together; according to EMWET 
    % documentation it must start at TE, then go over top to LE, then back 
    % along the lower side to TE; the first point of xylower is excluded as 
    % it is the same as the last of the reversed xyupper
    coords = [xyupper(length(xyupper):-1:1,:) % reverse
              xylower(2:length(xylower),:)]; % exclude first
    
    % write to file
    fid = fopen('withcomb135_EMWET.dat','wt');
    for i = 1:length(coords)
        fprintf(fid, '%g %g \n', coords(i,1), coords(i,2));
    end
    fclose(fid);
end


