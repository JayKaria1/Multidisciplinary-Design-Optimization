function MTOW = objective(X)
    Reference_Values;
    Initial_Geometry;
    
    fid = fopen('Design Vectors.txt','a');
    fprintf(fid, '%f\t', X);
    fprintf(fid, '\n');
    fclose(fid);

    global W_AminusW;    
    global CD0;
    
    global initial;

    W_fuel_guess = Ref.DFL;
    W_wing_guess = Ref.W_wing;
    
    [W_fuel,W_wing,counter] = MDA_coordinator(X,W_fuel_guess,W_wing_guess);
    
    MTOW_denorm = W_fuel+W_wing+W_AminusW;

    global couplings;  % Make a csv file that appends per iter.
    
    couplings.W_fuel = W_fuel;
    couplings.W_wing = W_wing;

    Init_MTOW = initial.W_fuel + initial.W_wing + W_AminusW;
    MTOW = MTOW_denorm/Init_MTOW;
    
    fid = fopen('iter_track.txt','a');
    fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\n', counter, MTOW, W_fuel, W_wing, W_AminusW, CD0);
    fclose(fid);



end