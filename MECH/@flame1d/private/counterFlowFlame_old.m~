function obj = counterFlowFlame_old(obj)
%obj = counterFlowFlame_old(obj)
% It needs a proper defined flame1d object, i.e you should set the
% compositions, the temperature, pressure etc. 
% Then you will enjoy a simulation of an axisymmetric stagnation-point 
% non-premixed flame. !!! OLD Version !!! Uses a reaction removed mechanism
% file for computiong the inital guess. Also, two internal gas object are
% created...
%
% "Life, loathe it or ignore it, you can't like it."
%		-- Marvin, the android,
%                    "Hitchhiker's Guide to the Galaxy"
% 
% (c) 2011 by Uwe Pruefert for VIRTUHCON
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 

obj.flameType = 'counterFlowFlame';
obj = obj.clearFractions;

% define the flame ...
% parameter values

if ~(isnumeric(obj.mdotLeft)&&isnumeric(obj.mdotRight))
    error('counterFlame:WrongInputFormat',...
                'The inlet flow velocities must be positive numbers.')
end
if ~(ischar(obj.compLeft)&&ischar(obj.compRight))
    obj.compLeft
    obj.compRight
    error('counterFlowFlame:MustBeChar',...
        'Wrong format of the composition in In- and/or Outlet given.')
end

if ~isempty(obj.mechanism)
    rxnmech = obj.mechanism;
else
    error('counterFlame:NoMechFile','Cannot found a mechanism file');
end

% identify the species given by the composition in the inlets
initialSpecies = [getSpeciesFromInitialFractions(obj.compLeft),...
    getSpeciesFromInitialFractions(obj.compRight)];

% the mechanism without any reactions should be nemaed by xxxwor.cti
rxnmechwor = [rxnmech(1:max(findstr('.',rxnmech))-1),'wor.cti'];

% create the do nothing mechanism

removeReactionsCTI(rxnmech,rxnmechwor);

% We assum, that the transpot with mix is "gas_1". Later we should look for
% a better solution, e.g. read the related line  in the cti file...
if isempty(obj.transport)
    transport  =  'gas_1'; 
else
    transport  = obj.transport;
end
    
           
% space up to now fix given
if isempty(obj.space)
    fprintf('The space property is empty, I will use a distance between the inlets of 0.02 meters.\n');
    fprintf('The grid is initialzed by 20 equidistant distributed points.\n');
    space = linspace(0,1,20);
    spaceStrech = 0.02;
    initial_grid = spaceStrech*space;  % m
else
    if length(obj.space)>=20&&length(obj.space)<=300  
        fprintf('Valid initial grid given.\n')      
        initial_grid = obj.space;
        space = linspace(0,1,length(initial_grid));
        spaceStrech = max(initial_grid);
    else
        fprintf('No valid initial grid given (too less or to manny points), but I have some data to create one.\n')
        if length(obj.space)>300
            initial_grid = linspace(obj.space(1),obj.space(end),300);
            space =  linspace(0,1,300);
            spaceStrech = max(initial_grid); 
        else
            initial_grid = linspace(obj.space(1),obj.space(end),20);
            space =  linspace(0,1,20);
            spaceStrech = max(initial_grid); 
        end           
    end
end
    
% hmmm, a good guess ???
tol_ss    = [1.0e-8 1.0e-12];       % [rtol atol] for steady-state problem
tol_ts    = [1.0e-6 1.0e-8];        % [rtol atol] for time stepping

loglevel  = 1;                      % amount of diagnostic output (0
                                    % to 5)
refine_grid = 0;                    % 1 to enable refinement, 0 to
                                    % disable ...	
                                    % It seems to be the best to swich it
                                    % off, because in view of a later
                                    % refinement, additional points can
                                    % *confuse* the solver 

				   
%%%%%%%%%%%%%%%% create the gas object %%%%%%%%%%%%%%%%%%%%%%%%
%
% This object will be used to evaluate all thermodynamic, kinetic,
% and transport properties
%
% 
% We need for the inital guess a non reacobj.initialTemperatureg gas  

gas0 = IdealGasMix(rxnmechwor,transport);
% equilibrate(gas0,'HP')
gas1 = IdealGasMix(rxnmech,transport);

% set its state to that of the  fuel (arbitrary)
set(gas1,'T', obj.initialTemperature, 'P', obj.initialPressure, 'X', obj.compRight);



%%%%%%%%%%%%%%%% create a flow object %%%%%%%%%%%%%%%%%%%%%%%
% 
% This flow is used by the guess and the real computations
%

f = AxisymmetricFlow(gas0,'flow');
f1 = AxisymmetricFlow(gas1,'flow');
set(f, 'P', obj.initialPressure, 'grid', initial_grid);
set(f, 'tol', tol_ss, 'tol-time', tol_ts);
set(f1, 'P', obj.initialPressure, 'grid', initial_grid);
set(f1, 'tol', tol_ss, 'tol-time', tol_ts);

%%%%%%%%%%%%%%% create the air inlet %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  The temperature, mass flux, and composition (relative molar) may be
%  specified.
%
% inletLeft = Inlet('air_inlet');
inletLeft = Inlet('inlet');
if isempty(obj.tInLeft)
    set(inletLeft, 'T', obj.initialTemperature, 'MassFlux', obj.mdotLeft, 'X', obj.compLeft);
else
    set(inletLeft, 'T', obj.tInLeft, 'MassFlux', obj.mdotLeft, 'X', obj.compLeft);
end
% sets the MOLE Fractions on the inlets


%%%%%%%%%%%%%% create the fuel inlet %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%
% inletRight = Inlet('fuel_inlet');
inletRight = Inlet('inlet');
if isempty(obj.tInRight)
    set(inletRight, 'T', obj.initialTemperature, 'MassFlux', obj.mdotRight, 'X', obj.compRight);
else
    set(inletRight, 'T', obj.tInRight , 'MassFlux', obj.mdotRight, 'X', obj.compRight);
end
% sets the MOLE Fractions on the inlets

%%%%%%%%%%%%% create the flame object  %%%%%%%%%%%%
%
% Once the component parts have been created, they can be assembled
% to create the flame objects.
% fl0 is for the guess
% fl1 for the final computation
%
fl0 = flame(gas0, inletLeft, f, inletRight);

% if obj.massFractions and obj.space fit together, we assume it as a
% initial guess
fits = ~isempty(obj.massFractions)&&(length(obj.massFractions(1,:))==length(obj.space))...
            &&(length(obj.temp)==length(obj.space));
if  fits
    % copy the mass fractions and the space to the flame
    yeq = obj.massFractions;   
    T1= obj.temp;
    zuV = [obj.space/spaceStrech;obj.u;obj.v];
    fprintf('\n\t\t\t\t#########################################################\n')
    fprintf('\t\t\t\t#     Inital guess taken  from values found in          #\n')
    fprintf('\t\t\t\t#     obj.temp, obj.massFractions, obj.u, and obj.v     #\n')
    fprintf('\t\t\t\t#########################################################\n\n\n') 
else
    %%%%%%%%%%%%%%%%  no initial guess, we try to compute one %%%%%%%%%
    % To compute the flux profile, we
    % solve without any reactions...
    %
    try
        solve(fl0, loglevel,refine_grid);         
        zuV = [z(f)/spaceStrech;solution(fl0,'flow','u');solution(fl0,'flow','V')];
        fprintf('\n\t\t\t\t#########################################################\n')
        fprintf('\t\t\t\t#                   Inital guess: compute               #\n')
        fprintf('\t\t\t\t#              solution of a pure flow problem          #\n')
        fprintf('\t\t\t\t#########################################################\n\n\n')       
    catch ME
        % Here we are really in trouble:
        fprintf('The transport equation should be solvable in every case so I\n')
        fprintf('stopp here all computations. Exeption message was:')        
        disp(ME.message)
        fprintf('Sorry for your inconvinience!')
        return   
    end
    
    obj = createSolution(obj,gas0,fl0,f);
%     obj.mixFractions = obj.getMixFractions;
    
   
    % Now we have the mixing of unburned gas
    % so we can compute the equilibrium in every point
    %
    for k = 1:length(initialSpecies)
        fractions(k,:) = solution(fl0,'flow',initialSpecies{k});
    end
        
    % z(f) returns the domain
    set(gas0,'Temperature',obj.initialTemperature,'Pressure',obj.initialPressure,'X', ...
            makeFractions(initialSpecies,fractions(:,k)));
%     yeq1(:,1) = massFractions(gas0);
%     T1(1) = temperature(gas0); 
    
    
    for k = 1:length(space)%      
        set(gas0,'Temperature',obj.initialTemperature,'Pressure',obj.initialPressure,'X', ...
            makeFractions(initialSpecies,fractions(:,k)));
        yeq0(:,k) = massFractions(gas0);
        T1(k) = temperature(gas0); 
        try
            equilibrate(gas0,'HP'); 
            yeq1(:,k) = massFractions(gas0);   
            T1(k) = temperature(gas0); 
        catch ME
            fprintf(['\nEquilibration failed in space point s(',num2str(k),')\n']);
            fprintf('Message was:\n\t\t\t\t');
            fprintf(ME.message)
            yeq1(:,k) =  yeq1(:,k-1);   
            T1(k) =  T1(k-1);
        end
        
    end
    %%%%% compute a approximation on the temperature in the flame
    s  = (T1-min(T1))/max(T1-min(T1));
    [~,indx] = max(s);
    leftOfFlame =obj.tInLeft*(space<space(indx)); 
    rightOfFlame =  obj.tInRight*(space>space(indx)); 
    T4  = s.*T1+ leftOfFlame+ rightOfFlame;
    % define the guess for the species when computed
    yeq = [yeq0(:,1),yeq1(:,2:end-1),yeq0(:,end)];
end


% the new flame object 
fl1 = flame(gas1, inletLeft, f1, inletRight);

%%%%%%%%%%%%%%%%%%% set the inital guess %%%%%%%%%
% eigther the computed one or the values in obj.massFractions
%  
%
setProfile(fl1, 2, 'T', [space;T4]);


% setProfile(fl1, 2, 'T', [space;T1]); %  fallback, when

% set flux and temperature profiles
setProfile(fl1, 2, {'u', 'V'},...
    zuV);
for n = 1:nSpecies(gas0)   
    nm = speciesName(gas0,n); 
    setProfile(fl1, 2, nm, [space;[massFraction(inletLeft,n),yeq(n,2:end-1),massFraction(inletRight,n)]])
end

fprintf('\n\t\t\t\t#########################################################\n')
fprintf('\t\t\t\t#                Mixture is now in equlibrium           #\n')
fprintf('\t\t\t\t#               by using equilibrate(gas,''HP'')          #\n')
fprintf('\t\t\t\t#########################################################\n\n\n')     

%%%%%%%%%%%% enable the energy equation %%%%%%%%%%%%%%%%%%%%%
%
%  The energy equation will now be solved to compute the
%  temperature profile. We also tighten the grid refinement
%  criteria to get an accurate final solution.
%
%


% enable the energie equation for flow f1

% setRefineCriteria(fl1, 1.4, 300.0, 0.1, 0.1);
try
    enableEnergy(f1);
    solve(fl1, loglevel,refine_grid);
    fprintf('\n\t\t\t\t#########################################################\n')
    fprintf('\t\t\t\t#          Solution computed with energy eq.            #\n')
    fprintf('\t\t\t\t#########################################################\n\n\n')  
catch ME
    % last try to solve the energy equation with the T1 profile
    fprintf('Solving the energie equation failed at first try, Exeption message was:') 
    disp(ME.message) 
    fprintf('Try now with the T1 temperature profile\n');
    try % ..a fall back
        setProfile(fl1, 2, 'T', [space;T1]);
        solve(fl1, loglevel,refine_grid)
    catch ME
        fprintf('\n\t\t\t\t#########################################################\n')
        fprintf('Solving the energie equation failed at second try, Exeption message was:') 
        disp(ME.message) 
        fprintf('Try now with solvig first problem without energy eqn.\n');
        fprintf('\n\t\t\t\t#########################################################\n\n\n')
        try
            setProfile(fl1, 2, 'T', [space;T1]);
            disableEnergy(f1)
            solve(fl1, loglevel,refine_grid);
            fprintf('\n\t\t\t\t#########################################################\n')
            fprintf('\t\t\t\t#        Solution computed without energy eq.           #\n')
            fprintf('\t\t\t\t#########################################################\n\n\n')  
        catch ME
            % Here we are  on the point where a failure can be caused be a bad
            % chemistry. We can fill the fild with zeros an teturn
            fprintf('\n\t\t\t\t#########################################################\n')
            fprintf('\t\t\t\t#  The integration failed, I stopp here all computations.\n#')
            fprintf('\t\t\t\t# Exeption message was:')  
            disp(ME.message)
            fprintf('\n\t\t\t\t#########################################################\n\n\n')
            obj.temp = zeros(size(z(f)));
            obj.space = z(f);
            obj = obj.setMassFractions(zeros(nSpecies(gas1),length(z(f))));
            obj = obj.setMoleFractions(zeros(nSpecies(gas1),length(z(f))));
            return    
        end
        enableEnergy(f1);
        solve(fl1, loglevel,refine_grid);
    end
end


%%%%%%%%%% make output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj = createSolution(obj,gas1,fl1,f1);

 

%%%%%%%%%% local hero fucntions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    function obj = createSolution(obj,gas1,fl1,f1)
        % we should first clear the old values, because we indexing the new ones by
        % a "row set call"
        mm = molarMasses(gas1);
        obj.temp = solution(fl1,'flow','T');
        obj.space = z(f1);
        sn = speciesNames(gas1);
        
        obj = obj.clearFractions;
        
        for kk = 1:nSpecies(gas1)   
             massFrac (kk,:) = solution(fl1,'flow',sn{kk});    
        end
        obj = obj.setMassFractions(massFrac);
        D = sparse(1:length(mm),1:length(mm),1./mm); 
        M = 1./(1./mm'*obj.massFractions);
         
        obj = obj.setMoleFractions((D*obj.massFractions)*sparse(1:length(M),1:length(M),M));
        obj.u = solution(fl1,'flow','u');
        obj.v = solution(fl1,'flow','V').*z(f1);
    end

end

