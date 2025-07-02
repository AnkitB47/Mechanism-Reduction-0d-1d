function  individualsListOut = evaluate(individualsListIn,refSolObj,objective,varargin)
%individualsListOut = evaluate(individualsListIn,referenceSolutionobject,evaluateFunction[,optional arguments])
global mechanismpath 
user = getenv('USER');

if isempty(mechanismpath)
    mechanismpath = ['/home/',user,'/mechanisms/'];
end


 
[status,result]  = system(['mkdir ',mechanismpath]);
if status~=0
    fprintf(1,['evaluate: A (possible) problem occured when trying mkdir. The message was:\n\t',...
        result,'\n']);    
end
% open the data file



fprintf(1,'^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^\n');
fprintf(1,'!This is evaluate, it evaluates (sic!) a generation of   species in the GA.!\n');
fprintf(1,'vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv\n');
% some input checking

if ~isa(individualsListIn,'individual')
    error('evaluate:InPutMustBeIndividual','The input must be a cell array. Each entry must be a individual object')
end


% a hiden folder for all mechanisms

% some definitions


if length(varargin)==1
    cases = 2;
    if ~isa(refSolObj,'referenceSolutionCube')
        error('evaluate:WrongClass','The reference solution object must be from class referenceSolutionCube.')
    end
    if isa(varargin{1},'scenario')
        scen = varargin{1};
        if isempty(scen.temperatures)
            scen = scen.setTemperatures(refSolObj.initialTemperature);
        end
        if isempty(scen.pressures)
            scen = scen.setPressures(refSolObj.initialPressure);
        end
        if isempty(scen.moleFractions)
            scen = scen.addMoleFractions(refSolObj.initialMoleFractions);
        end
    else
        error('evaluate:ArgumentNotSupported','The argument should be a scenario object.')
    end
elseif length(varargin)>1
    error('evaluate:TooManyInputs','Too many input arguments found.')
else
    if ~isa(refSolObj,'chemistry')
        error('evaluate:WrongClass','The reference solution object must be from class chemistry.')
    end
    cases = 1;
end
    

n = length(individualsListIn);


% the main loop over all individuals
for k = 1:n
    
    
    
    % switch the case "one szenario" or "multiszenario" evaluation
    switch cases
        case 1 % single scenario evaluation
            % if the value field is empty,
            % perform the evaluation defines by the function_handle
            if isempty(individualsListIn(k).value)
                % a name for the mechanism and create the mechanism file
                % for this individual
                fprintf(1,'\t\t\t\t------------------\n');
                fprintf(1,'\t\t\t\t! Remove species !\n');
                fprintf(1,'\t\t\t\t------------------\n');
                mechName = [mechanismpath,'mech',num2str(randi(100000000,1)),'.',refSolObj.mechanism(end-2:end)];
                removeSpecies(refSolObj.mechanism,mechName,individualsListIn(k).speciesNotInclude,...
                    individualsListIn(k).names(individualsListIn(k).speciesNotInclude));
                individualsListIn(k).mechanism = mechName;
                % prepare the chemistry object
                sol = chemistry(individualsListIn(k).mechanism);
                sol.initialMoleFractions = refSolObj.initialMoleFractions;
                sol.initialMassFractions = refSolObj.initialMassFractions;
                sol.compLeft = refSolObj.compLeft;
                sol.compRight = refSolObj.compRight;
                sol.mdotLeft = refSolObj.mdotLeft;
                sol.mdotRight = refSolObj.mdotRight;
                sol = sol.setSpace(refSolObj.space);
                sol.initialPressure = refSolObj.initialPressure;
                sol.initialTemperature = refSolObj.initialTemperature; 
                sol.simulationTime = refSolObj.simulationTime;
                % evaluation, solve the reactor and evaluate the objective
                sol = sol.setReactor(refSolObj.reactor);                           
                sol = feval(sol.reactor,sol);  
                
                %%%%%%%%%%%%%%%%%%%
%                 hold on
%                 if strcmp(refSolObj.reactor,'counterFlowFlame')
%                     plot(sol.space,sol.temp,'k');
%                     drawnow
%                 else
%                     plot(sol.times,sol.temp,'k')
%                     drawnow
%                 end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                individualsListIn(k).value = feval(objective,refSolObj,sol);
                [fid,message] = fopen([mechanismpath,'results.txt'],'a');
                if fid <= 2
                    fprintf(['evaluate. Something wrong with fopen. The message was: ',message,'\n'])
                end
                if individualsListIn(k).value~=Inf,
                    fprintf(fid,individualsListIn(k).mechanism);
                    fprintf(fid,'\tnspecies: %4d \tn of living species: %4d  ignition: %6.4e f(T,X) = %6.4e \n',...
                        length(getSpeciesFromMech(mechName)), sol.nLivingSpecies, sol.ignt, individualsListIn(k).value); 
                
                else
                    % eval is infty, we clear the mechfile, and don't write
                    % anything in the results.txt file
                    delete(mechName);
                end
                    fclose('all');
                    sol.clear;
               
            else
                fprintf(1,'Warning: individual has a value, skippig evaluation:\t');
                fprintf(1,[individualsListIn(k).mechanism,'\n']);
                [fid,message] = fopen([mechanismpath,'results.txt'],'a');
                if fid <= 2
                    fprintf(['evaluate. Something wrong with fopen. The message was: ',message,'\n'])
                end
                fprintf(fid,individualsListIn(k).mechanism);
                fprintf(fid,'\t\t\t\t\t');
                fprintf(fid,'f(T,X) = %6.4e',individualsListIn(k).value);
                fprintf(fid,'\n');
                fclose('all');
            end
            
            % some documentation
            
        case 2 % multi scenario evaluation
            mechName = [mechanismpath,'mech',num2str(randi(100000000,1)),'.',refSolObj.mechanism(end-2:end)];
            removeSpecies(refSolObj.valueCube(1,1,1).mechanism,mechName,individualsListIn(k).speciesNotInclude,...
                individualsListIn(k).names(individualsListIn(k).speciesNotInclude));
            % now, the individual get its mechanism
            individualsListIn(k).mechanism = mechName;
            nTemp = length(scen.temperatures);
            nPress= length(scen.pressures);
            nMoleFracs = length(scen.moleFractions);
            individualsListIn(k).value = 0;
            sol = chemistry(individualsListIn(k).mechanism);
            sol.simulationTime = refSolObj.valueCube(1,1,1).simulationTime;
            fprintf(fid,[mechName,': \n']);
            fprintf(fid,[mechName,'\n']);
            livingspecies = 0;
            for l = 1:nTemp
                sol.initialTemperature = scen.temperatures(l);
                for m = 1:nPress
                    sol.initialPressure = scen.pressures(m);
                    for r = 1:nMoleFracs
                        sol.initialMoleFractions = scen.moleFractions{r};
                        sol = feval(sol.reactor,sol);
                        % perform the evaluation defines by the
                        % function_handle evaluation
                        val= feval(objective,refSolObj.valueCube(l,m,r),sol); 
                        individualsListIn(k).value = ...
                            individualsListIn(k).value+val; 
                        fprintf(fid,'   nspecies: %4d nliving species: %4d ignition: %6.4e f(T,X) = %6.4e \n',...
                            [sol.nSpecies  sol.nLivingSpecies sol.ignt val]); 
                        fprintf(fid,'   nspecies: %4d nliving species: %4d ignition: %6.4e f(T,X) = %6.4e \n',...
                            [sol.nSpecies  sol.nLivingSpecies sol.ignt val]); 
                        if sol.nLivingSpecies >livingspecies
                            livingspecies = sol.nLivingSpecies;
                        end
                    end
                end
            end
            
            individualsListIn(k).value = individualsListIn(k).value*1/l/m/r;
            fprintf(fid,'\n\n------------------+++++++++++++++++++------------------\n\n');
            fprintf(fid,'   nspecies: %4d nliving species: %6.4e f(T,X) = %6.4e \n',...
                            [sol.nSpecies  livingspecies  individualsListIn(k).value]);  
            fprintf(fid,'\n\n------------------+++++++++++++++++++------------------\n\n');
            sol.clear;
        otherwise
    end
 
    
end

% creating the output
individualsListOut = individualsListIn;
[fid,message] = fopen([mechanismpath,'results.txt'],'a');
if fid <= 2
    fprintf(['evaluate. Something wrong with fopen. The message was: ',message,'\n'])
end
fprintf(fid,'\n\n------------------ xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ------------------\n\n');
fclose('all');
% local functions
   


end


