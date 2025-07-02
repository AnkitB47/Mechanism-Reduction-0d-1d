    classdef flame1d < chemistry    % 
        properties (SetAccess = private)                        
            mixFractions;
            u;
            v;
            transport;
        end
        
        properties (SetAccess = public)  
            space;
            compLeft;
            compRight;
            initialTemperature;
            initialPressure;
            mdotLeft;
            mdotRight;
            tInLeft;
            tInRight;
            flameType;        
        end
        
        
        methods(Access=public)
            % The CONSTRUCTOR
            function obj = flame1d(varargin)
                obj = obj@chemistry();
                switch nargin
                    case 0
                        % do nothing
                    case 2  
                        obj.setMechanism(varargin{1},varargin{2}); 
                        obj.transport = varargin{2};  
                    otherwise
                        ME = MException('WrongNumberArgin',...
                            ['Wrong number of arguments.',...
                            ' flame1d needs two arguments',...
                            'path to mechanism file and mixture model.']); 
                        rethrow(ME);
                end
            
                
            end
            
            
        
            % copy method for hardcopies        
            function obj2 = copy(obj1)
                if ~obj1.isvalid
                    ME = MException('batchReactor:copy',...
                        'Invalid or deleted object');
                    throw(ME)
                end
                obj2 = flame1d();            
                obj2.temp = obj1.temp;
                obj2.dens = obj1.dens;
                obj2.pres = obj1.pres;          
                obj2.mixFractions = obj1.mixFractions;
                obj2.u = obj1.u;
                obj2.v  = obj1.v;
                obj2.transport  = obj1.transport;
                obj2.space  = obj1.space;
                obj2.compLeft  = obj1.compLeft;
                obj2.compRight  = obj1.compRight;
                obj2.initialTemperature  = obj1.initialTemperature;
                obj2.initialPressure  = obj1.initialPressure;
                obj2.mdotLeft  = obj1.mdotLeft;
                obj2.mdotRight  = obj1.mdotRight;
                obj2.tInLeft  = obj1.tInLeft;
                obj2.tInRight  = obj1.tInRight;
                obj2.flameType  = obj1.flameType;                   
                obj2.mechanism  = obj1.mechanism;        
                obj2.gas  = obj1.gas;
                obj2.speciesNames  = obj1.speciesNames;
                obj2.massFractions  = obj1.massFractions;
                obj2.moleFractions   = obj1.moleFractions;            
            end
            
            function obj = setSpace(obj,space)
                % sets the private propertie chemistry.space
                % obj = obj.setSpace(space) where space can be a skalar, a
                % Vector with start and end value or a the discrete space
                % The space has to be in the form [0, s], hence setSpace will
                % transform every other intervall.   
                fprintf('\n');
                if ~isnumeric(space)
                    error('chemistry:spaceNumeric','The space must be numeric.')
                else
                    if isempty(space)
                        obj.space = [];
                    elseif min(size(space))>1
                        % matrix...
                        error('chemistry:spaceNumeric','The space must be a vector.')
                    elseif max(size(space))==1
                        % skalar: It will be interpreted as distance between
                        % the inlets...
                        obj.space = [0, space];
                        fprintf('Parameter s is given: Space is set to:\n')
                        disp(obj.space);
                    else
                        % The good case...
                        obj.space = space-space(1);
                        fprintf('A discrete space is given, will use it.\n')
                        if length(obj.space)>200,
                            fprintf('The inital discrete space has over 200 points. The can cause performance problems.\n')
                        else
                            %
                        end
                    end
                end
            end
            
            function solve(obj)
                switch obj.flameType
                    case 'counterFlowFlame'
                        obj.counterFlowFlame; 
                    case 'burnerFlame'
                        obj.burnerFlame;
                    otherwise
                end
            end
            
            % overwrite the methods from chemistry...
            function molefrac = getMoleFractionsByName(obj,Name)
                % obj.getMoleFractionsByName(Name)
                % returns the values of mole fractions over the time for species NAME
                % Name can be a string containing the species name or a cell
                % array of species names 
                if isa(Name,'cell')  
                    molefrac =  zeros(lenght(Name),...
                        length(obj.getMoleFractionsByName(Name{1})));         
                    for k = 1:length(Name)   
                        molefrac(k,:) = obj.getMoleFractionsByName(Name{k});
                    end
                elseif isa(Name,'char') 
                    molefrac = obj.moleFractions(obj.getSpeciesIndexByName(Name),:);
                    if isempty(molefrac),
                        warning('getMoleFractions:SpeciesNotFound',...
                            ['Species ',Name,' not found, is it in the mechanism?\n']); 
                        molefrac = zeros(size(length(obj.space)));
                    end
                else
                    error('getMoleFractionsByName:cellArray','Name must be a cell array or a string')
                end 
            end 
            
            function massfrac = getMassFractionsByName(obj,Name)
                % returns the values of mass fractions over the time for species NAME
                % Name can be a string containing the species name or a cell
                % array of species names 
                if isa(Name,'cell')                
                massfrac =  zeros(lenght(Name),...
                        length(obj.getMoleFractionsByName(Name{1}))); 
                    for k = 1:length(Name)
                        massfrac(k,:) = obj.getMassFractionsByName(Name{k});
                    end
                elseif isa(Name,'char') 
                    massfrac = obj.massFractions(obj.getSpeciesIndexByName(Name),:);
                    if isempty(massfrac),
                        warning('getMassFractionsByName:SpeciesNotFound',...
                            ['Species ',Name,' not found, is it in the mechanism?\n']); 
                        massfrac = zeros(size(length(obj.space)));
                    end
                else
                    error('getMassFractionsByName:cellArray','Name must be a cell array or a string')
                end 
            end                 
            
            function varargout = locateFlameFront(obj) 
                % gives back the location of the flame (front). 
                % [k,x] = obj.locateFlameFront  gives both x value and index
                % k = obj.locateFlameFront  gives only the index.
                % To determine x and k, we identify the flame(front) by the maximum
                % of the temperature.
                if isempty(obj.temp)
                    error('flame1d:locateFlameFront:emptyTemp',...
                        'Before computing the location of the flame you must compute the temperature.')
                else
                    [~,indx] = max(obj.temp);
                    switch nargout
                        case 0
                            fprintf(['\tFlame front located at x = ',num2str(obj.space(indx)),' , k = ',num2str(indx),'.\n'])
                        case 1
                            varargout{1} = indx;
                        case 2
                            varargout{1} = obj.space(indx);
                            varargout{2} = indx;
                        otherwise
                            error('flame1d:locateFlameFront:tooManyOutputs',...
                                'Too many output arguments, try >>help locateFlame.')
                    end                       
                end
            end
            
            function varargout = meanTemperature(obj)
                % t_mean =  obj.meanTemperature 
                % computes the mean temperature of a flame1d object
                % by t_mean = \int_x_0^x_1 T dx / \int_x_0^x_1 1 dx.
                if isempty(obj.temp) || isempty(obj.space)
                    error('flame1d:meanTemperature:emptyTemp',...
                        'Before computing the location of the flame you must compute the temperature.')
                else
                    h = obj.space(2:end)-obj.space(1:end-1);               
                    mT = 0.5/(max(obj.space)-min(obj.space))*h*(obj.temp(1:end-1)+obj.temp(2:end))';
                    switch nargout
                        case 0 
                            fprintf(['\t mean temperature = ',num2str(mT),'\n']);
                        case 1                        
                            varargout{1} = mT;                    
                        otherwise
                            error('flame1d:meanTemperature:tooManyOutputs',...
                                'Too many output arguments, try >>help meanTemperature.')
                    end                       
                end 
            end
            
        
            
            %%%%%%%%%% methods block ends here %%%%%%%%%%
        end
    end
