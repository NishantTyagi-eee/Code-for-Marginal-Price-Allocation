% This code is meant to compute the marginal prices of inertia and FFR
% services. The data has been pre-fed as per the Thesis. But the user may
% un-comment the input commannds for using this code with any general data.
clear all
clc

% STEP 1: Collect Data from User

%----------Step 1.1: Collecting information on Thermal Generators----------
% fprintf('------------------------------------------------------------------------------\n');
% fprintf('Firstly, we shall collect data on Thermal Generators of each TYPE.\n');
% fprintf('------------------------------------------------------------------------------\n');
% Therm_Gen_Min_Max=input('Please enter column vector containing the Min. and Max. PRODUCTION CAPACITY (in MW) for a Thermal Generators of each TYPE:\n'); 
% Therm_NoLoadCost=input('Please enter column vector containing NO-LOAD COSTS (in £) for a Thermal Generator of each TYPE:\n');
% Therm_MargCost=input('Please enter column vector containing MARGINAL COSTS (in £/MWh) for a Thermal Generator of each TYPE:\n');
% numThermGenPerType = input('Please enter column vector containing the NUMBER for a Thermal Generator of each TYPE:\n'); 
% Therm_H=input('Please enter column vector containing INERTIA CONSTANTS (in seconds) for a Thermal Generator of each TYPE:\n');
% T=input('Please enter column vector containing FR DELIVERY TIMES (in seconds) for a Thermal Generator of each TYPE:\n');
% T_Delay=input('Please enter column vector containing FR DELAY TIMES (in seconds) for a Thermal Generator of each TYPE:\n');
% HeadroomForFR=input('Please enter column vector containing the FRACTION (0 to 1) of AVAILABLE HEADROOM for FR from a Thermal Generator of each TYPE:\n');
% Headroom_Load_Slope=input('Please enter column vector containing the SLOPE (0 to 1) of how AVAILABLE HEADROOM for FR decreases with increasing LOAD from a Thermal Generator of each TYPE:\n');
%--------------------------------------------------------------------------

Therm_Gen_Min_Max=[0 1600; 0 1600]; % 2 types of Gen. one supplies (For BESS, Min_Max = 0 50 MW)
Therm_NoLoadCost=[750 700]'; % No load Cost_CCGT > No Load Cost_OCGT since CCGT have higher upfront and maintenance costs (For BESS, NL cost = 2lakh pounds ?? Check)
Therm_MargCost=[50 60]'; % Marg Cost_CCGT < Marg Cost_OCGT since CCGT are more efficient (For BESS marginal costs = DC market clearing price = 60 pound per MW)
numThermGenPerType=[8 1]'; % U can assume one battery
Therm_H=[4.5 3.7]'; % CCGT units have higher H since they are usually larger + operate at higher loads + include a steam turbine in addition to the gas turbine, which adds to the total rotating mass.
% OCGT units have lower H since are typically smaller as they used for peaking power or fast response situations. 
T=[20 2]';
T_Delay=[10 1.5]';
HeadroomForFR=[1 1]'; % What fraction of HEADROOM is AVAILABLE for FR 
Headroom_Load_Slope=[1 1]'; % How does that AVAIBLE HEADROOM for FR decrease with increasing LOAD
% Simultaneous clearing of energy and reserve has already been advocated by many ==> O.F. must be composed of costs of energy and costs of FR services
Therm_FR_Costs=0.001*ones(length(numThermGenPerType),1);
%('Please enter the Bids for providing FR service by Thermal Generators of each TYPE in £.\n','If the Bids are unknown please enter a small value like 0.01 for each'); 

%---------Step 1.2: Collecting information on Nuclear Generators-----------
% fprintf('------------------------------------------------------------------------------\n');
% fprintf('Secondly, we shall collect data on Nuclear Generators.\n');
% fprintf('------------------------------------------------------------------------------\n');
% Nuclear_Gen_Max=input('Please enter the PRODUCTION CAPACITY (in MW) for Nuclear Generator:\n'); 
% Nuclear_NoLoadCost=input('Please enter the NO-LOAD COSTS (in £) for Nuclear Generator:\n');
% Nuclear_MargCost=input('Please enter MARGINAL COSTS (in £/MWh) for Nuclear Generator:\n');
% Nuclear_H=input('Please enter the INERTIA CONSTANT (in seconds) for Nuclear Generator:\n');
% Nuclear_Part_Load=input('If the Nuclear Gen. is PART-LOADED, please enter the MW by which its output has been reduced\n');
% fprintf('The Nuclear Generator is not considered to provide FR due to technical limitations and economic considerations.\n');
% fprintf('Hence we are not collectiong data on FR Delivery Time for Nuclear Generator.\n');
%--------------------------------------------------------------------------
Nuclear_Gen_Max=1600;
Nuclear_MargCost=10;
Nuclear_H=5;
Nuclear_Part_Load=0;

%----------Step 1.3: Collecting information on Wind Plants-----------------
% fprintf('------------------------------------------------------------------------------\n');
% fprintf('Thirdly, we shall collect data on Wind Plants.\n');
% fprintf('------------------------------------------------------------------------------\n');
% Wind_Gen_Max=input('Please enter the PRODUCTION CAPACITY (in MW) for Wind Plants:\n'); 
% Wind_H=input('Please enter the INERTIA CONSTANT (in seconds) for Wind Plants:\n');
%--------------------------------------------------------------------------
Wind_Gen_Max=0;
Wind_H=0;

%----------Step 1.4: Collecting information on Grid Requirements-----------
% fprintf('------------------------------------------------------------------------------\n');
% fprintf('Lastly, we shall collect data on the Grid Requirements.\n');
% fprintf('------------------------------------------------------------------------------\n');
% Fnadir_max=input('Please enter the permissible FREQUENCY NADIR (in Hz) as per Grid Code:\n'); % Nadir Required by Grid Code in Hz
% RoCoF_max=input('Please enter the permissible RoCoF (in Hz/s) as per Grid Code:\n'); % RoCoF Required in Hz/s
% F_nominal=input('Please enter the Nominal Frequency as per Grid Code:\n'); % Nominal Frequency Required in Hz
% Load=input('Please enter the LOAD DEMAND (in MW) to be fulfilled by the Grid:\n');
%--------------------------------------------------------------------------
Fnadir_max=0.20; % Aim is to keep frequency within operational limits
RoCoF_max=1;
Load=9000;
F_nominal=50;
% The largest possible outage is set by the largest single plant i.e., nuclear plant
Max_Infeed_Loss=1600; 


% Input the parameters of the corresponding SFR model and triangular injection 
TR=3.3900;
zeta=0.4318;
w_n=0.4400;
alpha=1.5429;
w_r=0.3969;
phi=2.4365;
GSFR_Act_Integral=0.0930;
t_n=4.6097;
t1=1.7767;
t2=9.6928;
tri_start=1.5767;
tri_rise=0.5;
tri_peak=2.0767;
tri_fall=7.7710;
tri_end=9.8477;
tri_base=8.2710;
f_n=49.7063;
f_n_req=49.8140;

% STEP 2: Create "INTEGER" Optimisation Problem

%--------------Step 2.1: Declaring Decision Variables (DVs)----------------
%-----Step 2.1.1: Declaring Decision Variables (DVs) for Thermal Gen.------
Types_Therm_Gen=length(numThermGenPerType); % Counting HOW many TYPES of Therm Gen are there 
Therm_Yg = intvar(Types_Therm_Gen,1,'full'); % COLUMN DV Vector that represents HOW MANY Thermal Gen. of each TYPER are COMMITTED (integer variable).
Therm_Pg = sdpvar(Types_Therm_Gen,1,'full'); % COLUMN DV Vector that represents the TOTAL POWER produced by Thermal Gen. of each TYPE (real variable)
R = sdpvar(Types_Therm_Gen,1,'full'); % COLUMN DV Vector that represents the TOTAL FR provided by Thermal Gen. of each TYPE
assign(R,0.1);
% NOTE: 'R' is a DV because it is not a known variable before we solve the optimisation problem 
%-----Step 2.1.2: Declaring Decision Variables (DVs) for Nuclear Gen.------
Nuclear_Pg = sdpvar(1); %  Scalar DV that represents the POWER produced by the Nuclear Gen.
%-----Step 2.1.3: Declaring Decision Variables (DVs) for Wind Plants-------
Wind_Pg_Curtailed = sdpvar(1); %  Scalar DV that represents the amount of curtailed wind power.
% NOTE: Amount of Wind Power curtalied is for the optmisation problem to optimise
%--------Step 2.1.4: Declaring Decision Variables (DVs) for Grid-----------
Infeed_PLoss = sdpvar(1); % Scalar DV that represents the size of contingency; 
% We will assume the Infeed_Ploss = Max_Infeed_Loss = Nuclear_Gen_Max 
%------------------Step 2.2: Formulating Constraints-------------------------
%----------Step 2.2.1: Formulating Constraints on Thermal Gen.---------------
R_Max=HeadroomForFR.*Therm_Gen_Min_Max(:,2); % MAX FR that CAN be provided by a Therm. Gen. = FR_Limit of that Therm. Gen. (in fraction between 0 and 1) * MAX POSSIBLE OUTPUT of that Therm. Gen.
Therm_Pg_Max=Therm_Yg.*Therm_Gen_Min_Max(:,2); % MAX OUTPUT that CAN provided by a Therm. Gen. = Committed status (T_g = 0 or 1) * MAX POSSIBLE OUTPUT of that Therm. Gen.
Therm_Gen_Constraints = [0 <= Therm_Yg <= numThermGenPerType,... No of COMIITED Thermal Gen. of each TYPE must be between 0 and the TOTAL no of Thermal Gen. of that TYPE
                         Therm_Yg.*Therm_Gen_Min_Max(:,1) <= Therm_Pg <= Therm_Yg.*Therm_Gen_Min_Max(:,2),...% EQN. (40) -- GENERATION LIMITS
                         0 <= R <= Therm_Yg.*R_Max,... % EQN. (41) -- MAX FR that can be provided at any time = the MAX OUTPUT from COMMITTED Gen. at that time
                         R <= Headroom_Load_Slope.*(Therm_Pg_Max-Therm_Pg),...% EQN. (42) -- MAX FR that can be provided at any time = HEADROOM AVAILABLE at that time
                         0 <= R <= numThermGenPerType.*R_Max,...% ADDITIONAL CONSTRAINT on FR from Economic Dispatch Constraints Eqn. (30) (by LUIS BADESA) MAX FR that can be provided at any time = No of Therm. Gen. of a type * MAX FR from Therm. Gen. of that type
                         0 <= Therm_Pg <= numThermGenPerType.*Therm_Gen_Min_Max(:,2)]; % ADDITIONAL CONSTRAINT on Pg from Economic Dispatch Constraints Eqn. (29) (by LUIS BADESA) MAX Pg that can be provided at any time = No of Therm. Gen. of a type * MAX Pg from Therm. Gen. of that type
%---------Step 2.2.2: Formulating Constraints on Nuclear Gen.----------------
Nuclear_Gen_Constraints = (Nuclear_Gen_Max-Nuclear_Part_Load) <= Nuclear_Pg <= Nuclear_Gen_Max; % Nuclear Gen. Output lies between its DE-RATED CAPACITY and FULL CAPACITY 
%---------Step 2.2.3: Formulating Constraints on Wind Plants-----------------
Wind_Gen_Constraints = 0 <= Wind_Pg_Curtailed <= Wind_Gen_Max; % EQN. (43) -- RES Curtailment must be limited below RES production
%---------Step 2.2.4: Formulating Constraints on Infeed Loss-----------------
Infeed_Loss_Constraints = Infeed_PLoss == Nuclear_Pg;  % Infeed power loss = power dispatched from the nuclear unit.
% Additional Infeed_Loss_Constraint: (Max_Infeed_Loss-Nuclear_Part_Load) <= Infeed_PLoss <= Max_Infeed_Loss,...
%---------Step 2.2.5: Formulating Constraints on Power Balance--------------
Wind_Gen_Actual=Wind_Gen_Max-Wind_Pg_Curtailed;
P_Balance_Constraint = (sum(Therm_Pg) + Nuclear_Pg + Wind_Gen_Actual) == Load;
%-----------Step 2.2.6: Formulating RoCoF Constraint-------------------------
% To formulate RoCoF constraint, we need the TOTAL POST-FAULT inertia 
Total_H_PostFault = (((Therm_H.*Therm_Yg)'*Therm_Gen_Min_Max(:,2)) + (Wind_H*Wind_Gen_Actual))/F_nominal; % Eqn. (37) -- Therm_Inertia depends on Therm_Yg and we omit the Nuclear_H since it is disconnected
RoCoF_Constraint = Total_H_PostFault >= Infeed_PLoss/(2*RoCoF_max);% Eqn. (4) -- For RoCoF to be satisfied, the TOTAL POST-FAULT inertia must be => Infeed_PLoss/(2*RoCoF_Grid_limit)
%-----------Step 2.2.7: Formulating QSS Constraint---------------------------
QSS_Constraint = (Infeed_PLoss <= sum(R)); % Eqn. (5) -- For the frequency to stabilize in QSS, the TOTAL FR provided must be => Size of Disturbance
%----------Step 2.2.8: Formulating F_Nadir_Constraints-----------------------
%------------------ Formulation of F_Nadir_Constraints begins--------------   

% STEP (a): Finding the status of EACH FR in EACH of the time intervals
% (This step is quite different than the case of FR services WITHOUT activation delays.
% In that case, the status of FR services upto nth time interval was straightforward.
% The FR services upto n-1 were fully delivered, while FR services from n to S were ramping up.
% This facility was because the user entered the delivery periods of different FR services in ascending order.
% In case of FR services with activation delays, the t-axis is split in a rather complex manner.
% Hence, we need to make a deliberate attempt to find the status of FR services in each time interval.)
t_axis = unique([T_Delay, T]); % define time axis 
t_intervals = length(t_axis)-1; % No. of Time intervals is 1 less than the no of elements on time axis
%% We ignore the 1st time interval less than the LEAST T_Delay element since
%% nadir can't occur when no FR has been activated
S=length(T); % finding the number of FFR services
FR_status = zeros(S,t_intervals); % Defining a matrix 'FR_status' for the status of FR services 
% '0' = not yet started; '1' = yet ramping up; '2' = fully delivered;
% ROWS are the FR services; COLUMNS are the Time Intervals
for t=1:t_intervals
    t_begin=t_axis(t);
    t_end=t_axis(t+1);
       for k = 1:S
            if t_end < T_Delay(k)  % If the end point of (t)th t_interval lies BEFORE (k)th FR sevice BEGINS, its status is '0'
                FR_status(k,t) = 0;
            end
            if t_begin >= T_Delay(k) && t_end <= T(k) % If the start point of (t)th t_interval lies ON/AFTER (k)th FR service begins "AND" end point of (t)th t_interval lies ON/BEFORE (k)th FR sevice ENDS, its status is '1'
                FR_status(k,t) = 1;
            end
            if t_end > T(k) % If the end point of (t)th t_interval lies AFTER (k)th FR service ENDS, its status is '2'
                FR_status(k,t) = 2;
            end
       end
end
% STEP (b): Formulating the F_Nadir_Constraints & their Conditional Statements correponding to EACH of the time intervals  
    Compensated_H = sdpvar(1,t_intervals,'full'); % Declare a ROW DV VECTOR to store the 1st term on LHS (i.e., Y1 in Eqn 3.20 in PhD Thesis)
    FR_Ramping = sdpvar(1,t_intervals,'full'); % Declare a ROW DV VECTOR to store the 2nd term on LHS (i.e., Y2 in Eqn 3.20 in PhD Thesis)
    Compensated_PLoss = sdpvar(1,t_intervals,'full'); % Declare a ROW DV VECTOR to store the LHS Term in Eqn 3.20
    Accumulated_FR=sdpvar(1,t_intervals,'full'); % Declare a ROW DV VECTOR to store ACCUMULATED THERMAL FR delivered "UPTO" EACH of the time intervals
    % (The ACCUMULATED FR will be required for formulating the IF CONDITIONS associated with nadir contraints)
    for t=1:t_intervals
        t_end=t_axis(t+1);
        Compensated_H(t) = Total_H_PostFault; % Initialising Inertia Compensated by FR post-Fault
        FR_Ramping(t) = 0; % Initialising Acculumated FR at 0
        Compensated_PLoss(t) = Infeed_PLoss; % Initialising Power Loss Compensated by FR at the Post Fault power loss i.e., P_L
        Accumulated_FR(t)=0; % Initialising the Total FR accumulated upto EACH of the (t)th t_interval (needed for formulating IF CONDITIONS)
        for k = 1:S
            if FR_status(k,t) == 0 % kth FR service has not yet started delivering at (t)th t_interval
               Accumulated_FR(t) = Accumulated_FR(t)+0;
               Compensated_H(t) = Compensated_H(t) -0;
               Compensated_PLoss(t) = Infeed_PLoss - 0;
               FR_Ramping(t) = FR_Ramping(t)+0;
            end
            if FR_status(k,t) == 2 % kth FR service has "FULLY" delivered "UPTO" (t)th t_interval
                Accumulated_FR(t) = Accumulated_FR(t) + R(k); % Accumulating FR from FASTEST services which have been "FULLY DELIVERED" "UPTO" (t)th t_interval
                Compensated_H(t) = Compensated_H(t) - (R(k)*(T(k)+2*T_Delay(k))/(4*Fnadir_max)); % FASTEST services which have been "FULLY" delivered "UPTO" (t)th t_interval REDUCE overall system inertia requirement 
                Compensated_PLoss(t) = Compensated_PLoss(t) - R(k); % FASTEST services which have been "FULLY" delivered "UPTO" (t)th t_interval REDUCE power loss that remains to be compensated
            end
            if FR_status(k,t) == 1 % kth FR service is still ramping up at (t)th t_interval
                Accumulated_FR(t) = Accumulated_FR(t)+ R(k)*(t_end-T_Delay(k))/(T(k)-T_Delay(k)); % Accumulating FR from SLOWER services which are "STILL BEING" delivered at (t)th t_interval
                Compensated_H(t) = Compensated_H(t) + ((R(k)*((T_Delay(k))^2))/(T(k)*4*Fnadir_max)); % Contibution of services which are "STILL BEING" delivered at (t)th t_interval to Compensated_H
                Compensated_PLoss(t) = Compensated_PLoss(t) + (R(k)*T_Delay(k)/T(k)); % Contibution of services which are "STILL BEING" delivered at (t)th t_interval to Compensated_PLoss
                FR_Ramping(t) = FR_Ramping(t)+(R(k)/T(k)); % Accumulating FR from SLOWER services which are "STILL BEING" delivered at (t)th t_interval
            end
        end
    end

  % STEP (c): BOUNDING all the Terms in F_Nadir_Constraints because "implies function" requires this
  % "implies" function is implemented via Big-M approach ==> all its terms need to be EXPLICITY BOUNDED 
  % (Couldn't have known without LUIS BADESA)
 
    % BOUNDS on Compensated_H:
    LowerBound_Compensated_H = - numThermGenPerType'*(R_Max.*(T+2*T_Delay)/(4*Fnadir_max)) + numThermGenPerType'*(R_Max.*(T_Delay.^2)./(T*4*Fnadir_max));
    LowerBound_Compensated_H = LowerBound_Compensated_H*ones(1,t_intervals);
    UpperBound_Compensated_H = (((Therm_H.*numThermGenPerType)'*Therm_Gen_Min_Max(:,2)) + (Nuclear_H*Nuclear_Gen_Max) + (Wind_H*Wind_Gen_Max))/F_nominal + numThermGenPerType'*(R_Max.*(T_Delay.^2)./(T*4*Fnadir_max));
    UpperBound_Compensated_H = UpperBound_Compensated_H*ones(1,t_intervals);
    
    % BOUNDS on Compensated_PLoss:
    LowerBound_Compensated_PLoss = (Max_Infeed_Loss - numThermGenPerType'*R_Max + numThermGenPerType'*(R_Max.*T_Delay./T) );
    LowerBound_Compensated_PLoss = LowerBound_Compensated_PLoss*ones(1,t_intervals); 
    UpperBound_Compensated_PLoss = Max_Infeed_Loss + numThermGenPerType'*(R_Max.*T_Delay./T);
    UpperBound_Compensated_PLoss = UpperBound_Compensated_PLoss*ones(1,t_intervals); 
    
    % BOUNDS on FR_Ramping:
    LowerBound_FR_Ramping =0; % FR can never be negative
    UpperBound_FR_Ramping = numThermGenPerType'*(R_Max./T);
    UpperBound_FR_Ramping = UpperBound_FR_Ramping*ones(1,t_intervals);
    
    % BOUNDS on Accumulated_FR:
    LowerBound_Accumulated_FR =0; % FR can never be negative
    Upper_Bound_Accumulated_FR=numThermGenPerType'*R_Max;
    Upper_Bound_Accumulated_FR=Upper_Bound_Accumulated_FR*ones(1,t_intervals);
    
    % Collecting all the BOUNDS related to F_Nadir_constraints
    F_Nadir_Bounds=[LowerBound_Compensated_H <= Compensated_H <= UpperBound_Compensated_H,...
                    LowerBound_Compensated_PLoss <= Compensated_PLoss <= UpperBound_Compensated_PLoss,...
                    LowerBound_FR_Ramping <= FR_Ramping <= UpperBound_FR_Ramping,...
                    LowerBound_Accumulated_FR <= Accumulated_FR <= Upper_Bound_Accumulated_FR];

    % STEP (d): Implementing the "t_interval" F_Nadir_Constraints 
    F_Nadir_Constraints = []; % Initialising Matrix for F_Nadir_Constraints
    tol=0.0001; % Function "implies" doesn't work propely without a tolerance,
    for t=1:t_intervals
        if t==1 % For the F_Nadir_Constraint corresponding to FIRST interval (t=1), there is ONLY 1 condition to enforce 
                % It is that the Accumulated_FR "UPTO" 1st interval (Accumulated_FR(t=1)) should be > or = Infeed_PLoss minus a small tolerance (tol).
            F_Nadir_Constraints = [F_Nadir_Constraints,...
                implies(Accumulated_FR(t)>=Infeed_PLoss-tol,...
                [norm([Compensated_H(t)-FR_Ramping(t);2*sqrt(1/(4*Fnadir_max))*Compensated_PLoss(t)]);-Compensated_H(t);-FR_Ramping(t)]...
                <= [Compensated_H(t)+FR_Ramping(t);0;0])];
             % Why enforce non-negativity of Compensated_H and FR_Ramping? Since it is required for ROTATED SOC to be convex (LUIS BADESA)
        else % For the F_Nadir_Constraint corresponding to REMAINING time intervals (t>1), there are TWO conditions to enforce
             % CONDITION 1: Accumulated_FR "UPTO" (t-1)th i.e., PREVIOUS interval should be <= Infeed_PLoss + tolerance.
             % CONDITION 2: Accumulated_FR "UPTO" (t)th i.e., CURRENT interval should be >= Infeed_PLoss - tolerance. 
             F_Nadir_Constraints = [F_Nadir_Constraints,...
                implies([Accumulated_FR(t-1)<=Infeed_PLoss+tol,Accumulated_FR(t)>=Infeed_PLoss-tol],...
                [norm([Compensated_H(t)-FR_Ramping(t);2*sqrt(1/(4*Fnadir_max))*Compensated_PLoss(t)]);-Compensated_H(t);-FR_Ramping(t)]...
                <= [Compensated_H(t)+FR_Ramping(t);0;0])];
        end
    end
%------------------ Formulation of F_Nadir_Constraints ends----------------
% Collecting all the CONSTRAINTS together
Constraints = [Therm_Gen_Constraints,Nuclear_Gen_Constraints,Wind_Gen_Constraints,...
               Infeed_Loss_Constraints,P_Balance_Constraint,....
               QSS_Constraint,RoCoF_Constraint,F_Nadir_Bounds,F_Nadir_Constraints];
               
%Step 2.3: Declaring Objective Function
Objective = ((Therm_MargCost'*Therm_Pg) + (Therm_NoLoadCost'*Therm_Yg)) + (Nuclear_MargCost*Nuclear_Pg) + (Therm_FR_Costs'*R); %RES are assumed to have a zero marginal cost
%--------------------------------------------------------------------------

% STEP 3: Perform The INTEGER Optimization 
options = sdpsettings('solver','gurobi','gurobi.MIPGap',0.1e-2,'gurobi.QCPDual',1,'verbose',2); % Setting the Solver settings -- (Couldn't have known without LUIS BADESA)
sol = optimize(Constraints,Objective,options);
%--------------------------------------------------------------------------

% STEP 4: Processing and Displaying Results

%-----Step 4.1 Finding Time Interval in which F_Nadir_Constraint Occurs----
for t=1:t_intervals
    term1(t)=value(Compensated_H(t))*value(FR_Ramping(t));
    term2(t)=value(Compensated_PLoss(t))^2/(4*Fnadir_max);
    Nadir_Check(t)= term1(t) - term2(t);
    tolerance = 10;
    if ((Nadir_Check(t) > 0 || abs(Nadir_Check(t)) < tolerance) || (value(Accumulated_FR(t))>=value(Infeed_PLoss)))
    Nadir_Interval(t) = 1;
    else
    Nadir_Interval(t) = 0;
    end
   % Calculates a boolean value (1 for true, 0 for false) indicating whether F_Nadir occurs in t-th TIME INTERVAL
end
% if Nadir_Interval==0 && Nadir_Interval==1
%     fprintf('Nadir constraint is not satisfied');
% end

fprintf('------------------------------------------------------------------------------\n');
fprintf('The results of INTEGER OPTIMISATION run are as below:\n');
fprintf('The value of the objective function is:\n');
Optimum_Cost = value(Objective)*1e-3
fprintf('The Time Interval in which Frequency Nadir occurs is:\n');
Nadir_Time_Interval=find(Nadir_Interval)
fprintf('The Number of Units Committed for Thermal Generators of each TYPE is:\n');
Number_Committed_Per_Type=value(Therm_Yg)
fprintf('The Total Power (in MW) dispatched for Thermal Generators of each TYPE is:\n');
Thermal_Power_Dispatched_Per_Type = value(Therm_Pg)
fprintf('The Total FR (in MW) dispatched for Thermal Generators of each TYPE is:\n');
FR_Dispatched_Per_Type = value(R)
fprintf('The Total Power (in MW) dispatched for Wind Plants is:\n');
Wind_Power_Dispatched = Wind_Gen_Max-value(Wind_Pg_Curtailed)
fprintf('Next, performing the RELAXED OPTIMISATION run............:\n');
%--------------------------------------------------------------------------

% STEP 5: Create "RELAXED" Optimisation Problem
% Clearing all YALMIP variables but not the User input data
clear Therm_Yg Therm_Pg R Nuclear_Pg Wind_Pg_Curtailed Infeed_PLoss 
clear Accumulated_FR Compensated_H Compensated_PLoss FR_Ramping FR
clear Constraints Objective options sol

%--------------Step 5.1: Declaring Decision Variables (DVs)----------------
%-----Step 5.1.1: Declaring Decision Variables (DVs) for Thermal Gen.------
Types_Therm_Gen=length(numThermGenPerType); % Counting HOW many TYPES of Therm Gen are there 
Therm_Yg = sdpvar(Types_Therm_Gen,1,'full'); % COLUMN DV Vector that represents HOW MANY Thermal Gen. of each TYPER are COMMITTED (integer variable).
Therm_Pg = sdpvar(Types_Therm_Gen,1,'full'); % COLUMN DV Vector that represents the TOTAL POWER produced by Thermal Gen. of each TYPE (real variable)
R = sdpvar(Types_Therm_Gen,1,'full'); % COLUMN DV Vector that represents the TOTAL FR provided by Thermal Gen. of each TYPE
assign(R,0.1);
%-----Step 5.1.2: Declaring Decision Variables (DVs) for Nuclear Gen.------
Nuclear_Pg = sdpvar(1); %  Scalar DV that represents the POWER produced by the Nuclear Gen.
%--------------------------------------------------------------------------
%-----Step 5.1.3: Declaring Decision Variables (DVs) for Wind Plants-------
Wind_Pg_Curtailed = sdpvar(1); %  Scalar DV that represents the amount of curtailed wind power.
%--------------------------------------------------------------------------
%--------Step 5.1.4: Declaring Decision Variables (DVs) for Grid-----------
Infeed_PLoss = sdpvar(1); % Scalar DV that represents the size of contingency; 
% We will assume the Infeed_Ploss = Max_Infeed_Loss = Nuclear_Gen_Max 
%--------------------------------------------------------------------------

%------------------Step 5.2: Formulating Constraints-------------------------
%----------Step 5.2.1: Formulating Constraints on Thermal Gen.---------------
R_Max=HeadroomForFR.*Therm_Gen_Min_Max(:,2); % MAX FR that CAN be provided by a Therm. Gen. = FR_Limit of that Therm. Gen. (in fraction between 0 and 1) * MAX POSSIBLE OUTPUT of that Therm. Gen.
Therm_Pg_Max=Therm_Yg.*Therm_Gen_Min_Max(:,2); % MAX OUTPUT that CAN provided by a Therm. Gen. = Committed status (T_g = 0 or 1) * MAX POSSIBLE OUTPUT of that Therm. Gen.
Therm_Gen_Constraints = [0 <= Therm_Yg <= numThermGenPerType,... No of COMIITED Thermal Gen. of each TYPE must be between 0 and the TOTAL no of Thermal Gen. of that TYPE
                         Therm_Yg.*Therm_Gen_Min_Max(:,1) <= Therm_Pg <= Therm_Yg.*Therm_Gen_Min_Max(:,2),...% EQN. (40) -- GENERATION LIMITS
                         0 <= R <= Therm_Yg.*R_Max,... % EQN. (41) -- MAX FR that can be provided at any time = the MAX OUTPUT from COMMITTED Gen. at that time
                         R <= Headroom_Load_Slope.*(Therm_Pg_Max-Therm_Pg),...% EQN. (42) -- MAX FR that can be provided at any time = HEADROOM AVAILABLE at that time
                         0 <= R <= numThermGenPerType.*R_Max,...% ADDITIONAL CONSTRAINT on FR from Economic Dispatch Constraints EQN. (30) (by LUIS BADESA) MAX FR that can be provided at any time = No of Therm. Gen. of a type * MAX FR from Therm. Gen. of that type
                         0 <= Therm_Pg <= numThermGenPerType.*Therm_Gen_Min_Max(:,2)]; % ADDITIONAL CONSTRAINT on Pg from Economic Dispatch Constraints EQN. (29) (by LUIS BADESA) MAX Pg that can be provided at any time = No of Therm. Gen. of a type * MAX Pg from Therm. Gen. of that type
%--------------------------------------------------------------------------
%---------Step 5.2.2: Formulating Constraints on Nuclear Gen.----------------
Nuclear_Gen_Constraints = (Nuclear_Gen_Max-Nuclear_Part_Load) <= Nuclear_Pg <= Nuclear_Gen_Max; % Nuclear Gen. Output lies between its DE-RATED CAPACITY and FULL CAPACITY 
%--------------------------------------------------------------------------
%---------Step 5.2.3: Formulating Constraints on Wind Plants-----------------
Wind_Gen_Constraints = 0 <= Wind_Pg_Curtailed <= Wind_Gen_Max; % EQN. (43) -- RES Curtailment must be limited below RES production
%--------------------------------------------------------------------------
%---------Step 5.2.4: Formulating Constraints on Infeed Loss-----------------
Infeed_Loss_Constraints = Infeed_PLoss == Nuclear_Pg;  % Infeed power loss = power dispatched from the nuclear unit.
% Additional Infeed_Loss_Constraint: (Max_Infeed_Loss-Nuclear_Part_Load) <= Infeed_PLoss <= Max_Infeed_Loss,...
%--------------------------------------------------------------------------
%---------Step 5.2.5: Formulating Constraints on Power Balance--------------
Wind_Gen_Actual=Wind_Gen_Max-Wind_Pg_Curtailed;
P_Balance_Constraint = (sum(Therm_Pg) + Nuclear_Pg + Wind_Gen_Actual) == Load;
%--------------------------------------------------------------------------
%-----------Step 5.2.6: Formulating RoCoF Constraint-------------------------
% To formulate RoCoF constraint, we need the TOTAL POST-FAULT inertia 
Total_H_PostFault = (((Therm_H.*Therm_Yg)'*Therm_Gen_Min_Max(:,2)) + (Wind_H*Wind_Gen_Actual))/F_nominal; % EQN. (37) -- Therm_Inertia depends on Therm_Yg and we omit the Nuclear_H since it is disconnected
RoCoF_Constraint = Total_H_PostFault >= Infeed_PLoss/(2*RoCoF_max);% EQN. (4) -- For RoCoF to be satisfied, the TOTAL POST-FAULT inertia must be => Infeed_PLoss/(2*RoCoF_Grid_limit)
%-----------Step 5.2.7: Formulating QSS Constraint---------------------------
QSS_Constraint = sum(R) >= Infeed_PLoss; % EQN. (5) -- For the frequency to stabilize in QSS, the TOTAL FR provided must be => Size of Disturbance
%----------Step 5.2.8: Formulating F_Nadir_Constraints-----------------------
%------------------ Formulation of F_Nadir_Constraints begins--------------   
% STEP (a): Finding the status of EACH FR in EACH of the time intervals
t_axis = unique([T_Delay, T]); % define time axis 
t_intervals = length(t_axis)-1; % No. of Time intervals is 1 less than the no of elements on time axis
% We ignore the 1st time interval less than the LEAST T_Delay element since
% nadir can't occur when no FR has been activated
S=length(T); % finding the number of FFR services
FR_status = zeros(S,t_intervals); % Defining a matrix 'FR_status' for the status of FR services 
% '0' = not yet started; '1' = yet ramping up; '2' = fully delivered;
% ROWS are the FR services; COLUMNS are the Time Intervals
for t=1:t_intervals
    t_begin=t_axis(t);
    t_end=t_axis(t+1);
       for k = 1:S
            if t_end < T_Delay(k)  % If the end point of (t)th t_interval lies BEFORE (k)th FR sevice ENDS, its status is '0'
                FR_status(k,t) = 0;
            end
            if t_begin >= T_Delay(k) && t_end <= T(k) % If the start point of (t)th t_interval lies ON/AFTER (k)th FR service begins "AND" end point of (t)th t_interval lies ON/BEFORE (k)th FR sevice ENDS, its status is '1'
                FR_status(k,t) = 1;
            end
            if t_end > T(k) % If the end point of (t)th t_interval lies AFTER (k)th FR service ENDS, its status is '2'
                FR_status(k,t) = 2;
            end
       end
end
% STEP (b): Formulating the F_Nadir_Constraints & their Conditional Statements correponding to EACH of the time intervals  
    Compensated_H = sdpvar(1,t_intervals,'full'); % Declare a ROW DV VECTOR to store the 1st term on LHS (i.e., Y1 in Eqn 3.20 in PhD Thesis)
    FR_Ramping = sdpvar(1,t_intervals,'full'); % Declare a ROW DV VECTOR to store the 2nd term on LHS (i.e., Y2 in Eqn 3.20 in PhD Thesis)
    Compensated_PLoss = sdpvar(1,t_intervals,'full'); % Declare a ROW DV VECTOR to store the LHS Term in Eqn 3.20
    Accumulated_FR=sdpvar(1,t_intervals,'full'); % Declare a ROW DV VECTOR to store ACCUMULATED THERMAL FR delivered "UPTO" EACH of the time intervals
    % (The ACCUMULATED FR will be required for formulating the IF CONDITIONS associated with nadir contraints)
    for t=1:t_intervals
        t_end=t_axis(t+1);
        Compensated_H(t) = Total_H_PostFault; % Initialising Inertia Compensated by FR post-Fault
        FR_Ramping(t) = 0; % Initialising Acculumated FR at 0
        Compensated_PLoss(t) = Infeed_PLoss; % Initialising Power Loss Compensated by FR at the Post Fault power loss i.e., P_L
        Accumulated_FR(t)=0; % Initialising the Total FR accumulated upto EACH of the (t)th t_interval (needed for formulating IF CONDITIONS)
        for k = 1:S
            if FR_status(k,t) == 0 % kth FR service has not yet started delivering at (t)th t_interval
               % DO Nothing
            end
            if FR_status(k,t) == 2 % kth FR service has "FULLY" delivered "UPTO" (t)th t_interval
                Accumulated_FR(t)=Accumulated_FR(t) + R(k); % Accumulating FR from FASTEST services which have been "FULLY DELIVERED" "UPTO" (t)th t_interval
                Compensated_H(t) = Compensated_H(t) - (R(k)*(T(k)+2*T_Delay(k))/(4*Fnadir_max)); % FASTEST services which have been "FULLY" delivered "UPTO" (t)th t_interval REDUCE overall system inertia requirement 
                Compensated_PLoss(t) = Compensated_H(t) - R(k); % FASTEST services which have been "FULLY" delivered "UPTO" (t)th t_interval REDUCE power loss that remains to be compensated
            end
            if FR_status(k,t) == 1 % kth FR service is still ramping up at (t)th t_interval
                Accumulated_FR(t)=Accumulated_FR(t)+(R(k)*(t_end-T_Delay(k))/(T(k)-T_Delay(k))); % Accumulating FR from SLOWER services which are "STILL BEING" delivered at (t)th t_interval
                Compensated_H(t) = Compensated_H(t) + ((R(k)*((T_Delay(k))^2))/(T(k)*4*Fnadir_max)); % Contibution of services which are "STILL BEING" delivered at (t)th t_interval to Compensated_H
                Compensated_PLoss(t) = Compensated_H(t) + (R(k)*T_Delay(k)/T(k)); % Contibution of services which are "STILL BEING" delivered at (t)th t_interval to Compensated_PLoss
                FR_Ramping(t)= FR_Ramping(t)+(R(k)/T(k)); % Accumulating FR from SLOWER services which are "STILL BEING" delivered at (t)th t_interval
            end
        end
    end

  % STEP (c): BOUNDING all the Terms in F_Nadir_Constraints because "implies function" requires this
  % "implies" function is implemented via Big-M approach ==> all its terms need to be EXPLICITY BOUNDED 
  % (Couldn't have known without LUIS BADESA)
 
    % BOUNDS on Compensated_H:
    LowerBound_Compensated_H = - numThermGenPerType'*(R_Max.*(T+2*T_Delay)/(4*Fnadir_max)) + numThermGenPerType'*(R_Max.*(T_Delay.^2)./(T*4*Fnadir_max));
    LowerBound_Compensated_H = LowerBound_Compensated_H*ones(1,t_intervals);
    UpperBound_Compensated_H = (((Therm_H.*numThermGenPerType)'*Therm_Gen_Min_Max(:,2)) + (Nuclear_H*Nuclear_Gen_Max) + (Wind_H*Wind_Gen_Max))/F_nominal + numThermGenPerType'*(R_Max.*(T_Delay.^2)./(T*4*Fnadir_max));
    UpperBound_Compensated_H = UpperBound_Compensated_H*ones(1,t_intervals);
    
    % BOUNDS on Compensated_PLoss:
    LowerBound_Compensated_PLoss = (Max_Infeed_Loss - numThermGenPerType'*R_Max + numThermGenPerType'*(R_Max.*T_Delay./T) );
    LowerBound_Compensated_PLoss = LowerBound_Compensated_PLoss*ones(1,t_intervals); 
    UpperBound_Compensated_PLoss = Max_Infeed_Loss + numThermGenPerType'*(R_Max.*T_Delay./T);
    UpperBound_Compensated_PLoss = UpperBound_Compensated_PLoss*ones(1,t_intervals); 

    % BOUNDS on FR_Ramping:
    LowerBound_FR_Ramping =0; % FR can never be negative
    UpperBound_FR_Ramping = numThermGenPerType'*(R_Max./T);
    UpperBound_FR_Ramping = UpperBound_FR_Ramping*ones(1,t_intervals);
    
    % BOUNDS on Accumulated_FR:
    LowerBound_Accumulated_FR =0; % FR can never be negative
    Upper_Bound_Accumulated_FR=numThermGenPerType'*R_Max;
    Upper_Bound_Accumulated_FR=Upper_Bound_Accumulated_FR*ones(1,t_intervals);
    
    % Collecting all the BOUNDS related to F_Nadir_constraints
    F_Nadir_Bounds=[LowerBound_Compensated_H <= Compensated_H <= UpperBound_Compensated_H,...
                    LowerBound_Compensated_PLoss <= Compensated_PLoss <= UpperBound_Compensated_PLoss,...
                    LowerBound_FR_Ramping <= FR_Ramping <= UpperBound_FR_Ramping,...
                    LowerBound_Accumulated_FR <= Accumulated_FR <= Upper_Bound_Accumulated_FR];


    % STEP (d): Implementing the "t_interval" F_Nadir_Constraints 
    F_Nadir_Constraints = []; % Initialising Matrix for F_Nadir_Constraints
    tol=0.00001; % Function "implies" doesn't work propely without a tolerance,
    t=Nadir_Time_Interval; % We need to calculate ONLY 1 F_Nadir_Constraint corresponding to the Time Interval in which F_Nadir occurs 
        if t==1 % For the F_Nadir_Constraint corresponding to FIRST interval (t=1), there is ONLY 1 condition to enforce 
                % It is that the Accumulated_FR "UPTO" 1st interval (Accumulated_FR(t=1)) should be > or = Infeed_PLoss minus a small tolerance (tol).
                % We use "CONE" instead of "NORM" since we need to extract % VECTOR DUAL MULTIPLIERS -- (LUIS BADESA)
                % Cone has the syntax cone(x,y) <==> |x|<=y.
            F_Nadir_Constraints = [F_Nadir_Constraints,...
                                  (Accumulated_FR(t)>=Infeed_PLoss-tol):'limits',...
                                  (cone([Compensated_H(t)-FR_Ramping(t);2*sqrt(1/(4*Fnadir_max))*Compensated_PLoss(t)],Compensated_H(t)+FR_Ramping(t))):'nadir',... 
                                  Compensated_H(t)>=0,FR_Ramping(t)>=0];
            
            % Why enforce non-negativity of Compensated_H and FR_Ramping? Since it is required for ROTATED SOC to be convex           
        else % For the F_Nadir_Constraint corresponding to REMAINING time intervals (t>1), there are TWO conditions to enforce
             % CONDITION 1: Accumulated_FR "UPTO" (t-1)th i.e., PREVIOUS interval should be <= Infeed_PLoss + tolerance.
             % CONDITION 2: Accumulated_FR "UPTO" (t)th i.e., CURRENT interval should be >= Infeed_PLoss - tolerance. 
            F_Nadir_Constraints = [F_Nadir_Constraints,...
                                (Accumulated_FR(t-1)<=Infeed_PLoss+tol):'limits',...
                                (Accumulated_FR(t)>=Infeed_PLoss-tol):'limits',...
                                (cone([Compensated_H(t)-FR_Ramping(t);2*sqrt(1/(4*Fnadir_max))*Compensated_PLoss(t)],Compensated_H(t)+FR_Ramping(t))):'nadir',...
                                Compensated_H(t)>=0,FR_Ramping(t)>=0];
        end

        %------------------ Formulation of F_Nadir_Constraints ends----------------

Constraints = [Therm_Gen_Constraints,Nuclear_Gen_Constraints,Wind_Gen_Constraints,...
               Infeed_Loss_Constraints,P_Balance_Constraint,....
               QSS_Constraint,RoCoF_Constraint,F_Nadir_Bounds,F_Nadir_Constraints];
               
%Step 5.3: Declaring Objective Function
Objective = ((Therm_MargCost'*Therm_Pg) + (Therm_NoLoadCost'*Therm_Yg)) + (Nuclear_MargCost*Nuclear_Pg) + (Therm_FR_Costs'*R); %RES are assumed to have a zero marginal cost
%--------------------------------------------------------------------------

% STEP 6: Perform The RELAXED Optimization 
options = sdpsettings('solver','gurobi','gurobi.MIPGap',0.1e-2,'gurobi.QCPDual',1,'verbose',0); % Setting the Solver settings:
sol = optimize(Constraints,Objective,options);
%--------------------------------------------------------------------------


% STEP 7: Processing and Displaying Results
%-----------------Step 7.1 Extracting the DUAL VARIABLES-------------------
RoCoF_Dual_Variable=dual(RoCoF_Constraint);
QSS_Dual_Variable=dual(QSS_Constraint);
Nadir_Vector_Dual_Variables=dual(Constraints('nadir'));
mu=Nadir_Vector_Dual_Variables(1); % "mu" must be positive because of the DUAL FEASIBILITY Condition EQN. 22 (LUIS BADESA)
lambda1= -Nadir_Vector_Dual_Variables(2); % We inverse signs of λ1 and λ2 to satisfy COMPLEMENTARY SLACKNESS condition 
lambda2= -Nadir_Vector_Dual_Variables(3); % We inverse signs of λ1 and λ2 to satisfy COMPLEMENTARY SLACKNESS condition

%---Step 7.2 Calculating prices if COMPLEMENTARY SLACKNESS condition for Nadir Constraints is satisfied----
fprintf('Checking if the COMPLEMENTARY SLACKNESS condition is satisfied for the SOC Nadir Constraints.....\n');
t=Nadir_Time_Interval;
Comp_Slack_Check= lambda1*(value(Compensated_H(t))-value(FR_Ramping(t))) + lambda2*(2*sqrt(1/(4*Fnadir_max))*value(Compensated_PLoss(t)))...
                  - mu*(value(Compensated_H(t))+value(FR_Ramping(t))); % -- EQN. 21 ==> If this expression = 0, condition is satisfied 
tolerance = 5e-1; % tolerance to check if the value of Comp_Slack_Check is close enough to zero
if abs(Comp_Slack_Check) < tolerance % Checking the COMPLEMENTARY SLACKNESS condition for Nadir Constraints
    fprintf('The COMPLEMENTARY SLACKNESS condition for the SOC Nadir Constraints is satisfied.\n');
    fprintf('------------------------------------------------------------------------------\n');
    fprintf('The results of RELAXED OPTIMISATION run are as below:\n');
%   fprintf('The values of the DUAL VARIABLES are:\n');
%   RoCoF_Dual_Variable
%   QSS_Dual_Variable
%   Nadir_Vector_Dual_Variables
    fprintf('The Price of Inertia is:\n');
    Price_of_Inertia=(2*RoCoF_Dual_Variable)+(mu-lambda1)/F_nominal % Price in £/MWs
    fprintf('FASTER FR Services are those which are fully delivered by f_nadir\n');
    fprintf('SLOWER FR Services are those which are still ramping up by f_nadir\n');
    fprintf('\n');
    % Calculating and displaying the prices of FASTER FR Services
    for t=Nadir_Time_Interval % Looping through time intervals UPTO the time interval in which f_nadir occurs
        for k=1:S % Looping through the status of ALL the frequency services for one particualr time interval
            if FR_status(k,t)==2 % the (k)th FR service has been fully delivered by f_nadir
                fprintf('The price of FASTER FR Service Number %d (fully delivered by f_nadir) is:\n',k);
                Price_Faster_FR_Service = QSS_Dual_Variable - ((mu-lambda1)*(T(k)+(2*T_Delay(k)))/(4*Fnadir_max)) + lambda2/sqrt(Fnadir_max)
            end
        end
    end

     % Calculating and displaying the prices of SLOWER FR Services
    for t=Nadir_Time_Interval % Looping through time intervals UPTO the time interval in which f_nadir occurs
        for k=1:S % Looping through the status of ALL the frequency services for one particualr time interval
            if FR_status(k,t)==1 || FR_status(k,t)==0 % the (k)th FR service is still ramping up by f_nadir
                fprintf('The price of SLOWER FR Service Number %d (still ramping up by f_nadir) is:\n',k);
                Price_Slower_FR_Service = QSS_Dual_Variable + ((mu-lambda1)*((T_Delay(k))^2)/(4*T(k)*Fnadir_max)) + ((mu+lambda1)/T(k)) - ((lambda2*T_Delay(k))/(T(k)*sqrt(Fnadir_max)))
            end
        end
    end
else
    fprintf('The COMPLEMENTARY SLACKNESS condition for the SOC Nadir Constraints is NOT satisfied.\n');
end

