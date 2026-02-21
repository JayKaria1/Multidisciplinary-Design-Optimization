%This file is purely for testing function files. Please commment other files's tests and add the file that you would like to test.
%Please do not remove any lines of pre-written codes from this file

%Calling the Initial Geometry function file and reference values
% Initial_Geometry;
Reference_Values;
% disp("Initial run done")

%Test for Q3D solver: Aero function file has been called.
%Result : Successful 
%Values:CD = 0.0068, CL = 0.3752, LD = 54.8480, CD_AminusW = 0.0166,time = 0.6s approx on 22/3/2025
%Recent run values: LD = 25.3439, CL = 0.3891, CD = 0.0154 , CD_AminusW = 0.0090, time = 12.9 s approx on 26/3/2025
% Geom.Root_chord = (Area+(Inboard_span^2 + (Inboard_span*(Geom.Halfspan-Inboard_span)))*tand(QC_Sweep-TE_sweep_adjustment))/(2*(Inboard_span+(Geom.Halfspan-Inboard_span)*(1+TPR)/2));

% [LD,CL,CD] = Aero(Ref.Span,Ref.Area,Ref.TR,Ref.Sweep_25,Ref.UCST,Ref.LCST,Ref.MTOW,Ref.DFL,0);
% CD_AminusW = (CL/16)-CD

% Test for Q3D Solver: Loads function file
%Result: Successful
%Values: 
% [Loads_out] = Loads(Initial_span,Initial_area,Initial_TPR,Initial_QCsweep,UCST,LCST,Ref.MTOW, Ref.MaxSpeed)
 
% Test for Performance function file
%Result: Successful. Returns value for W_fuel
%Values:  W_fuel = 79757 [kg]
% [Perfomance_Test] = Performance(Ref.M_cruise, Ref.Alt, Ref.MTOW, LD)

%Test for Structures
%Result: Successful. Returns value for W_wing
%Values: W_wing = 19408 [kg]
% W_wing= Structures(Ref.DFL,Ref.Span,Ref.Area,Ref.TR,Ref.Sweep_25,Ref.UCST,Ref.LCST,20000,0.2,0.7,0) %W wing initial guess require

%Test for MDA coordinator (convergence loop)
%Result Successful. Gives the outputs as needed and responds to changes in
%design vector  X = [UCST1 -> UCST6, LCST1->LCST6, Span, Area, TPR, QC_Sweep, M_cruise, H_cruise, Spar_front, Spar_back]
X = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];

[W_fuel,W_wing] = Initial_run(X,Ref.DFL,Ref.W_wing)

% [W_fuel,W_wing] = MDA_coordinator(X,Ref.DFL,Ref.W_wing)
% % 
% global couplings
% couplings.W_fuel = W_fuel;
% couplings.W_wing = W_wing;
% 
% %Test for Constraints File. 
% %Result: Successful. 
% [c,ceq] = Constraints(X)