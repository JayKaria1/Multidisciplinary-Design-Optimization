clear all;
%File to run the optimization

%Initializing text files for output
fid1 = fopen('iter_track.txt', "w");
fclose(fid1);

fid5 = fopen('CST_show.txt',"w");
fclose(fid5);
fid = fopen('iter_track.txt','a');
fprintf(fid, 'Counter\tMTOW\tW_fuel\tW_wing\tW_a_minus_W_wing\tCD0\n');
fclose(fid);

fid3 = fopen('tolcheck.txt', "w");
fclose(fid3);

fid2 = fopen('tolcheck.txt', 'a');
fprintf(fid2, 'c1\tc2\tc3\n');
fclose(fid2);

fid4 = fopen('Design Vectors.txt','w');
fclose(fid4);



X0 = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
ub = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1];
lb = [-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1];
%Options for optimizer

options = optimoptions(@fmincon,...
        'Algorithm','sqp',...
        'Display','iter-detailed',...
        'ConstraintTolerance',1e-3,...
        'DiffMinChange',5e-2,...
        'DiffMaxChange',5e-1,...
        'PlotFcn',{'optimplotfval','optimplotx','optimplotfirstorderopt','optimplotconstrviolation'});

%Following lines of code run the optimizer. First refrence and initial
%geometry is loaded, then initial run is done before starting the
%optimization

Reference_Values;

tic
[W_fuel_init,W_wing_init] = Initial_run(X0,Ref.DFL,Ref.W_wing);
toc

disp('Starting Optimization')
tic
[X_opt,MTOW_min_denormalised,exitflag,output] = fmincon(@(X) objective(X),X0,[],[],[],[],lb,ub,@(X) Constraints(X),options);
toc

save('Optimized_X',"X_opt");




