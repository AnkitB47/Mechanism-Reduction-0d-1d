function pureMixing(obj)
%  pureMixing(obj)
% It needs a proper defined chemistry object, i.e you should set the
% compositions, the temperature, pressure etc. 
% then it computes the mixing of a non-primixed counter flow flame
% by convection and difffusion.
%  
% 
% (c) 2011 by Uwe Pruefert for VIRTUHCON
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj.clearFractions;

% define the flame ...
% parameter values

if ~(isnumeric(obj.mdotLeft)&&isnumeric(obj.mdotRight))
    error('counterFlame:WrongInputFormat',...
                'The inlet flow velocities must be positive numbers.')
end
if ~(ischar(obj.compLeft)&&ischar(obj.compRight))
    error('counterFlowFlame:MustBeChar',...
        'Wrong format of the composition in In- and/or Outlet given.')
end
if isempty(obj.gas)||~isIdealGas(obj.gas)
    try % to setup the gas here
        if isempty(obj.mechanism)            
            error('counterFlame:NoMechFile','Cannot found a mechanism file');
        end
        obj.setGas(IdealGasMix(obj.mechanism,'gas_1'));      
    catch ME
        error('counterFlowFlame:NoGas',...
            ME.message);
    end
% else
    % obj.gas is a idealGas so we can use it
end



 
           
% space up to now fix given
if isempty(obj.space)
    fprintf('The space property is empty, Set it to standard value: Use a distance between the inlets of 0.02 meters.\n');
    fprintf('The grid is initialzed by 20 equidistant distributed points.\n');
    space = linspace(0,1,20);
    spaceStrech = 0.02;
    initial_grid = spaceStrech*space;  % m
else
    if length(obj.space)>=20&&length(obj.space)<=300  
        fprintf('Valid initial grid given.\n')      
        initial_grid = obj.space;
%         space = linspace(0,1,length(initial_grid));
%         spaceStrech = max(initial_grid);
    else
        fprintf('No valid initial grid given (too less or to manny points), but I have some data to create one.\n')
        if length(obj.space)>300
            initial_grid = linspace(obj.space(1),obj.space(end),300);
%             space =  linspace(0,1,300);
%             spaceStrech = max(initial_grid); 
        else
            initial_grid = linspace(obj.space(1),obj.space(end),20);
%             space =  linspace(0,1,20);
%             spaceStrech = max(initial_grid); 
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

				   



% create a non-reaction gas by setting the reaction multiplier to
% zero
setMultiplier(obj.gas,1:nReactions(obj.gas),0);

% set its state to that of the  fuel (arbitrary)
set(obj.gas,'T', obj.initialTemperature, 'P', obj.initialPressure, 'X', obj.compRight);



%%%%%%%%%%%%%%%% create a flow object %%%%%%%%%%%%%%%%%%%%%%%
% 
% This flow is used by the guess and the real computations
%

f = AxisymmetricFlow(obj.gas,'flow');

set(f, 'P', obj.initialPressure, 'grid', initial_grid);
set(f, 'tol', tol_ss, 'tol-time', tol_ts);



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
fl0 = flame(obj.gas, inletLeft, f, inletRight);

% if obj.massFractions and obj.space fit together, we assume it as a
% initial guess
try
    solve(fl0, loglevel,refine_grid);
    fprintf('\n\t\t\t\t#########################################################\n')
    fprintf('\t\t\t\t#              Inital guess: computed solution          #\n')
    fprintf('\t\t\t\t#              of a pure flow problem                   #\n')
    fprintf('\t\t\t\t#########################################################\n\n\n')       
catch ME
    % Here we are really in trouble:
    fprintf('The transport equation should be solvable in every case so I\n')
    fprintf('stopp here all computations.')  
    fprintf('Sorry for your inconvinience!')
    error('pureMixing:Failed',['Exeption message was:',ME.message])  
end
createSolution(obj,fl0,f);
% now set all reaction multipliers back to one
 


%%%%%%%%%% make output %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 


%%%%%%%%%% local hero fucntions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    function createSolution(obj,fl1,f1)
        % we should first clear the old values, because we indexing the new ones by
        % a "row set call"
        mm = molecularWeights(obj.gas);
        obj.temp = solution(fl1,'flow','T');
        obj.space = z(f1);
        sn = speciesNames(obj.gas);
        
        obj.clearFractions;
        
        for kk = 1:nSpecies(obj.gas)   
            massFrac(kk,:) = solution(fl1,'flow',sn{kk});    
        end
        D = sparse(1:length(mm),1:length(mm),1./mm); 
        M = 1./(1./mm'*massFrac);
        obj.setMassFractions(massFrac);
        obj.setMoleFractions((D*massFrac)*sparse(1:length(M),1:length(M),M));
        obj.u = solution(fl1,'flow','u');
        obj.v = solution(fl1,'flow','V').*z(f1);
    end

end

