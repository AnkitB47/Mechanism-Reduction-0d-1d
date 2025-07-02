classdef chemistry < handle
    % Container class for chemical reactions. 
    % It is an  AbstracCLass and is the mother class for all
    % chemical reactor and flame classes.
    % The code uses the last version of MATLAB OO. 
    % Be sure to use MATLAB >R2011b.
    %
    % Version 2.0.1. (C) 2013 by U. Pruefert 
    % Methods:
    %
    %
    % getLivingSpeciesNumber    double  number of species with mass/mole
    %                                   fraction > eps
    % mole2mass                 double  mole fraction to mass fractions
    % setMechanism              self    sets mechanism for gas object
    %             
    % clearFractions            self    clear obj.mass/moleFractions
    % getMassFractionsByName    double  gets the mass fraction for selected species   
    % nActiveReactions          double  gives the number of active reaction
    %                                   in multiplier lambda
    % setMoleFractions          self    sets the mole fractions to value
    % copy                      chemistry copies the object (abstract method)
    % getMoleFractionsByName    double  gives the mole fractions for
    %                                   selected species
    % nInactiveReactions        double  number of in-active reaction in
    %                                   multiplier lambda
    % setPressure               self    set the pressure to value
    % delete                            deletes the object, not the variable                   
    % getReactionsMultiplier    double  gives multiplier for reactions
    % nReactions                double  gives number of reactions              
    % setReactionsMultiplier    self    sets reaction multiplier lambda
    % getDeadSpeciesIndex       double  index of dead species in species vector
    % getSpeciesIndexByName     double  index of given species in species vector 
    % nSpecies                  double  number of species            
    % setTemperature            self    sets temperature 
    % getDeadSpeciesNames       gives   dead species names
    % isvalid                   logical from handle
    % removeSpecies                     removes selected species from
    %                                   mechanism file
    % getDeadSpeciesNumber      double  gives number of dead species
    % lookForSpeciesPattern     indx    gives species that fits the name
    %                                   pattern
    % setDensity                self    set density to value
    % getLivingSpeciesIndex     double  gives the index of living species
    % mainSpecies               double  gives the index of main species
    % setGas                    self    sets the gas property
    % getLivingSpeciesNames     cell    gives names strings as cell array
    % mass2mole                 double  comuted mole fractions from mass
    %                                   fractions
    % setMassFractions          self    sets the mass fraction to value
    % 
    % Static methods:
    % 
    % getSpeciesFromMech        cell    gives the species names of all
    %                                   species in mechanism
    % isInList                  logical checks if species in species list
    % mixFractions              fractions mix mole/mass fractions
    % indexCtiFile              double  computes structure of mechanism        
    % isInListI                 logical checks if species in species list
    % vector2CanteraString      fractions computes cantera fractions string
    %   
    % For details see manual.
    
    
    
    
    
    
    
    properties (SetAccess = protected)
        mechanism;        
        gas;
        speciesNames;
        massFractions;
        moleFractions; 
        temp;
        dens;
        pres;   
        vol;
    end   
    
    properties (Constant,Access = protected)
        % define some errors objects
         wrongInputFileFormatException = MException('chemistry:WrongInputFileFormat',...
             'The input must be a cti, inp or xml File.');
         wrongInputClassException = MException('chemistry:WrongInputClassFormat',...
             'The input class is wrong. Check it.');                       
         wrongNumberInputException = MException('chemistry:WrongNumberInputs',...
             'Wrong number of input arguments. Check it.')
         wrongNumberOutputException = MException('chemistry:WrongNumberOutputs',...
             'Wrong number of output arguments. Check it.')
         wrongInputArgumentException = MException('chemistry:WrongInputArgument',...
             'The input argument is wrong. Check it.'); 
         noStringException = MException('chemistry:noStringException',...
             'The argument must be a string');       
         noFileException = MException('chemistry:noFileExceptionFound',...
             'Cannot open the mechanism file.');
         emptyFieldException = MException('chemistry:emptyField',...
             'The argument is empty.');
    end
    
    methods(Access = public)
        function obj = chemistry(varargin) 
            % class constructor
            % call with zero, one or two argumens
            % obj = chemistry()
            % obj = chemistry('/path/to/mechanism/file')
            % obj = chemistry('/path/to/mechanism/file','mixModel')
            % (c) 2013 by Uwe Prüfert
            switch nargin  
                case 0 
                    % allowes empty call
                    return
                case 1                    
                    try
                        obj.mechanism = varargin{1};
                        switch obj.mechanism(end-2:end)
                            case {'cti' 'inp' 'xml'}
                                try
                                    obj.gas = Solution(obj.mechanism); 
                                catch ME
                                    ME.throwAsCaller
                                end
                                obj.speciesNames = speciesNames(obj.gas); ...
                                    %#ok<CPROP> this ist the speciesName method of Solution
                            otherwise
                                throwAsCaller(obj.wrongInputFileFormatException);
                        end
                    catch ME                        
                        obj.gas = [];
                        obj.mechanism = [];
                        obj.speciesNames = [];
                        cleanup;
                        rethrow(ME)                        
                    end                    
                otherwise
                    throwAsCaller(obj.wrongInputFileFormatException)
            end
        end   
    
         
        % all the nice methods, but only useful if you have a child class of
        % chemistry
        function setGas(obj,SolutionObject)
            % obj.setGas(CanteraSolutionObject)
            % links the Solution-Object with obj.gas
            % handle with care!
            if isIdealGas(SolutionObject)
                obj.gas = SolutionObject;
            else
                error('chemistry:setGas:NoIdealGas',...
                    'Setting the gas failed: No Ideal Gas given.')
            end
        end
        
        % A set method. Since obj.mechanism is private, one need such
        % method to change the mechanism in an EXISTING Object.
        function setMechanism(obj,mechfile,mix)  
        %  obj.setMechanism(mechfile)        
        %  obj.setMechanism(mechfile,mixModell)    
        % A set method. Since obj.mechanism is private, one need such
        % method to change the mechanism in an EXISTING Object. 
            switch nargin
                case 1
                    % use the empty argument to remove all data 
                    % leave the gas living, is it a good idea???
                    obj.mechanism = [];
                    obj.speciesNames = [];                    
                case 2                 
                    if ~ischar(mechfile)
                        throwAsCaller(obj.noStringException)                          
                    end                   
                    try                           
                        obj.gas = Solution(mechfile);                                          
                    catch ME
                        warning('chemistry:setMechanism:CreatingGasFailed',...
                            'Creating the gas property failed. Error message was:')
                        fprintf(ME.message);
                        obj.gas = [];
                        obj.mechanism = [];
                        obj.speciesNames = [];
                        cleanup;
                        return
                    end
                    obj.mechanism = mechfile;
                    obj.speciesNames = speciesNames(obj.gas);%#ok<CPROP>
                    if isempty(obj.speciesNames)
                        error('chemistry:setMechanism:EMptySpeciesList',...
                            ['The species list is empty, prob.',...
                            ' an error while creating the gas object.'])     
                    end
                case 3
                    if ~ischar(mechfile)
                        throwAsCaller(obj.noStringException) 
                    end                   
                    try
                        obj.gas = Solution(mechfile,mix); 
                        obj.mechanism = mechfile;
                        obj.speciesNames = speciesNames(obj.gas);%#ok<CPROP>                       
                    catch ME
                        warning('chemistry:setMechanism:CreatingGasFailed',...
                            'Creating the gas property failed. Error message was:')
                        fprintf(ME.message);
                        obj.gas = [];
                        obj.mechanism = [];
                        obj.speciesNames = [];
                        cleanup;
                        return
                    end
                    obj.mechanism = mechfile;
                    obj.speciesNames = speciesNames(obj.gas);%#ok<CPROP>
                    if isempty(obj.speciesNames)
                        error('chemistry:setMechanism:EMptySpeciesList',...
                            ['The species list is empty, prob.',...
                            ' an error while creating the gas object.'])
                    end
                otherwise
                    error('setMechanism:TooManyArguments',...
                        'Too many arguments when calling setMechanism.')
            end
        end
       
        function indx = getSpeciesIndexByName(obj,varargin)
            % indx = obj.getSpeciesIndexByName(string)
            % indx = obj.getSpeciesIndexByName(species1,spec2,...,speciesN)
            % gives the index of a species in obj.speciesNames
            % Name can be a string containing the species name or a cell
            % array of species names. 
            
            if nargin==1
                obj.wrongNumberInputException.throwAsCaller;
            elseif nargin==2
                Name = varargin{1};
            else
                Name = varargin(:);
            end
            
            if isa(Name,'char')
                % single "char" input action
                indx = [];
                for k = 1:length(obj.speciesNames)
                    if strcmp(Name,obj.speciesNames{k})
                        indx = k;
                    end
                end
                if isempty(indx),
                    ME = MException('getSpeciesIndex:SpeciesNotFound',...
                        ['Species ',Name,' not found.']);
                    ME.throwAsCaller;
                elseif length(indx)>1
                    ME = MException('getSpeciesIndex:DublicateSpecies',...
                        ['Species ',Name,' found more then one times'...
                        'in the mechanism, something wrong?\n']);
                    ME.throwAsCaller;
                end
            elseif isa(Name,'cell')
                % cell array action
                indx = zeros(1,length(Name));
                for k = 1:length(Name)
                    % Recursion, baby!
                    indx(k) = obj.getSpeciesIndexByName(Name{k});
                end
            end
        end        
        
        % two function from gas
        function n = nSpecies(obj)
            % n = obj.nSpecies
            % the number of species in the object
            try
                n = nSpecies(obj.gas);
            catch ME
                ME.throwAsCaller;
            end
        end       
         
        function nreactions = nReactions(obj)   
            % n = obj.nReactions
            % returns the number of reactions in the mechanism
            nreactions = nReactions(obj.gas);
        end
        
        function setMoleFractions(obj,fractions)
            % obj.setMoleFractions(fractions)
            % must be a CANTERA fractions string like
            % 'CH4:1,H2O:0.1,CO2:0.05,N2:8' 
            switch nargin
                case 1
                    obj.moleFractions = [];
                case 2
                    obj.moleFractions = fractions;
                otherwise
                    throwAsCaller(chemistry.wrongNumberInputException)
            end
        end
        
        function setMassFractions(obj,fractions)
            % obj.setMassFractions(fractions)
            % must be a CANTERA fractions string like
            % 'CH4:1,H2O:0.1,CO2:0.05,N2:8' 
            switch nargin
                case 1
                    obj.massFractions = [];
                case 2
                    obj.massFractions = fractions;
                otherwise
                    throwAsCaller(chemistry.wrongNumberInputException)
            end
        end
        
        function obj = setTemperature(obj,temp)
            obj.temp = temp;
        end
        
        function obj = setDensity(obj,dens)
            obj.dens = dens;
        end
        function obj = setPressure(obj,pres)
            obj.pres = pres;
        end
        
        function clearFractions(obj)
            % obj.clearFractions 
            % overwtrites all Fractions
            obj.massFractions = [];
            obj.moleFractions = [];
        end
        
        function molefrac = getMoleFractionsByName(obj,varargin)
            % molefrac = obj.getMoleFractionsByName(Name)
            % obj.getMoleFractionsByName(string)
            % obj.getMoleFractionsByName(species1,spec2,...,speciesN)
            % returns the values of mole fractions over the time for
            % species NAME
            % Name can be a string containing the species name or a cell
            % array of species names 
            molefrac = [];            
            if nargin==1
                obj.wrongNumberInputException.throwAsCaller;
            elseif nargin==2
                Name = varargin{1};
            else
                Name = varargin(:);
            end
            if isa(Name,'cell')  
                molefrac   =  zeros(length(Name),length(obj.times)) ;                
                for k = 1:length(Name)                      
                    molefrac(k,:) = obj.getMoleFractionsByName(Name{k});
                end
            elseif isa(Name,'char')                     
                if  ~isempty(obj.moleFractions) 
                    molefrac = obj.moleFractions...
                        (obj.getSpeciesIndexByName(Name),:);
                end
            else
                throwAsCaller(obj.wrongInputClassException);
            end            
        end 
        
        function massfrac = getMassFractionsByName(obj,varargin)
            % massfrac = obj.getMassFractionsByName(Name)
            % returns the values of mass fractions over the time for 
            % species NAME
            % Name can be a string containing the species name or a cell
            % array of species names 
            massfrac = [];            
            if nargin==1
                obj.wrongNumberInputException.throwAsCaller;
            elseif nargin==2
                Name = varargin{1};
            else
                Name = varargin(:);
            end            
            if isa(Name,'cell')  
                massfrac   =  zeros(length(Name),length(obj.times)) ;                
                for k = 1:length(Name)                      
                    massfrac(k,:) = obj.getMassFractionsByName(Name{k});
                end
            elseif isa(Name,'char')
                if  ~isempty(obj.massFractions) 
                    massfrac = obj.massFractions...
                        (obj.getSpeciesIndexByName(Name),:);
                end                    
            else
                throwAsCaller(obj.wrongInputClassException);
            end             
        end         
        
        % some methods to hande dead species
        function indx = getDeadSpeciesIndex(obj,epsilon)
            % dead is a species iff massfraction = 0 for all times
            % or epsilon if an argument is given to getDeadSpeciesIndex
            if nargin == 2 && ~(epsilon <1 &&  epsilon >=0)
                ME = MException ('chemistry:InvalidData',...
                        'epsilon must be between 0 and 1.');
                    ME.throwAsCaller;                 
            elseif nargin == 1                
                epsilon = eps;
            end   
            if isempty(obj.moleFractions)&& ~isempty(obj.massFractions)
                mass = true;
            elseif ~isempty(obj.moleFractions)&& isempty(obj.massFractions)
                mass = false;
            elseif ~(isempty(obj.moleFractions)||isempty(obj.massFractions))
                mass = true;
            else
                ME = MException('chemistry:NoData',...
                         'chemistry.moleFractions is empty or not well formed.');         
                     ME.throwAsCaller
            end
                
            
            
            indx = [];
            for k = 1:length(obj.speciesNames)
                if mass
                    if length(obj.massFractions(:,1))~=obj.nSpecies
                         ME = MException('chemistry:NoData',...
                             [ 'Number of values in fractions vector '...
                             'is not equal the number of species.']);         
                         ME.throwAsCaller
                    end
                    if max(obj.massFractions(k,:))<=epsilon; 
                        indx = [indx k]; 
                    end
                else
                    if length(obj.moleFractions(:,1))~=obj.nSpecies
                         ME = MException('chemistry:NoData',...
                             [ 'Number of values in fractions vector '...
                             'is not equal the number of species.']);         
                         ME.throwAsCaller
                    end
                    if max(obj.moleFractions(k,:))<=epsilon; 
                        indx = [indx k]; 
                    end
                end
            end
        end        
        
        function n = getDeadSpeciesNumber(obj,epsilon)
            % the number of dead species
            % n = obj.getDeadSpeciesNumber()
            % n = obj.getDeadSpeciesNumber(epsilon)
            % "dead" is a species iff massfraction = 0 for all times
            % or leser that epsilon if an argument is given to getDeadSpeciesIndex
            if nargin == 2 && ~(epsilon <1 &&  epsilon >=0)
                ME = MException ('chemistry:InvalidData',...
                        'epsilon must be between 0 and 1.');
                    ME.throwAsCaller;                 
            elseif nargin == 2
                % do nothing (we want to use the argument epsilon)
            else
                epsilon = eps;
            end                    
            n = length(getDeadSpeciesIndex(obj,epsilon));            
        end        
        
        function dead = getDeadSpeciesNames(obj,epsilon)
            % species = obj.getDeadSpeciesNames()
            % species = obj.getDeadSpeciesNames(epsilon)
            % "dead" is a species iff massfraction = 0 for all times
            % or leser that epsilon if an argument is given to getDeadSpeciesIndex
            % species contains the names of dead species 
            if nargin == 2 && ~(epsilon <1 &&  epsilon >=0)
                ME = MException ('chemistry:InvalidData',...
                        'epsilon must be between 0 and 1.');
                    ME.throwAsCaller;                 
            elseif nargin == 2
                % do nothing (we want to use the argument epsilon)
            else
                epsilon = eps;
            end     
            indx = getDeadSpeciesIndex(obj,epsilon);
            dead = cell(1,length(indx));
            for k = 1:length(indx)
                dead{k} = obj.speciesNames{indx(k)};
            end
        end
        
        function indx = getLivingSpeciesIndex(obj,epsilon)
            % Index vector of the living species.
            % indx = obj.getLivingSpeciesIndex()
            % indx = obj.getLivingSpeciesIndex(epsilon)
            % "Dead" is a  species with  Mole fractions lesser than 1e-6.
            %  All other are living. Method 
            %  gives back the related index vector.
            if nargin == 2 && ~(epsilon <1 &&  epsilon >=0)
                ME = MException ('chemistry:InvalidData',...
                        'epsilon must be between 0 and 1.');
                    ME.throwAsCaller;                 
            elseif nargin == 2
                % do nothing (we want to use the argument epsilon)
            else
                epsilon = eps;
            end      
            dead = getDeadSpeciesIndex(obj,epsilon);             
            indx = 1:nSpecies(obj);            
            indx(dead) = [];
        end
        
         % some methods to hande living species
        function alive = getLivingSpeciesNames(obj,epsilon) 
            % the names of living species 
            % indx = obj.getLivingSpeciesNames()
            % indx = obj.getLivingSpeciesNames(epsilon)
            % "Dead" is a  species with  Mole fractions lesser than 1e-6.
            %  All other are living. Method 
            if nargin == 2 && ~(epsilon <1 &&  epsilon >=0)
                ME = MException ('chemistry:InvalidData',...
                        'epsilon must be between 0 and 1.');
                    ME.throwAsCaller;                 
            elseif nargin == 2
                % do nothing (we want to use the argument epsilon)
            else
                epsilon = eps;
            end     
            indx = getLivingSpeciesIndex(obj,epsilon);            
            alive = cell(1,length(indx));
            for k = 1:length(indx)
                alive{k} = obj.speciesNames{indx(k)};
            end
        end        
       
        function n = getLivingSpeciesNumber(obj,epsilon)
            % the number of living species             
            % n = obj.getLivingSpeciesNumber()
            % n = obj.getLivingSpeciesNumber(epsilon)
            % "dead" is a species iff massfraction = 0 for all times
            % or leser than epsilon if given.
            if nargin == 2 && ~(epsilon <1 &&  epsilon >=0)
                ME = MException ('chemistry:InvalidData',...
                        'epsilon must be between 0 and 1.');
                    ME.throwAsCaller;                 
            elseif nargin == 2
                % do nothing (we want to use the argument epsilon)
            else
                epsilon = eps;
            end     
            n = length(getLivingSpeciesIndex(obj,epsilon));             
        end      
        
        function mp = getReactionsMultiplier(obj)
            % gives back the vector of reaction multipliers
            mp = multiplier(obj.gas,1:obj.nReactions);
        end
        % some methods that handle reactions
        function n = nActiveReactions(obj)
            % the number of active reactions, i.e. where multiplier is
            % not equal to zero
            n = sum(obj.getReactionsMultiplier~=0);
        end
           
        function setReactionsMultiplier(obj,reactions)
            % method that switch reactions on/off
            % Example: setReactionMultiplier(obj,[1 0 1 ...])   
            % Vector must be a logical or numerical vector of length
            % nReactions.
            % Note, that we manipulate the handle_object
            % obj.gas so we did not need an output argument.
            reactions = int32(reactions); 
            indx = find(~reactions);
            setMultiplier(obj.gas,indx,0);
            fprintf(['Remove ',num2str(length(indx)),...
                ' reactions...\n']);            
            indx = find(reactions);
            setMultiplier(obj.gas,indx,1);
            fprintf(['leave ',num2str(length(indx)),...
                ' reactions in mechanism.\n']);% 
            
        end
        
        function n = nInactiveReactions(obj)
            % the number of inactive reactions, i.e. where multiplier is
            % equal to zero
            n = sum(obj.getReactionsMultiplier==0);
        end
        
        % species handling
        function [varargout] = lookForSpeciesPattern(obj,varargin)
            % lookForSpeciesPattern(pattern)
            % Gives all species back where patter is found
            % Example
            %>> [~,sp] = sol.lookForSpeciesPattern('CH2')
            %
            %sp = 
            %
            %    'CH2'    'CH2(S)'    'CH2O'    'CH2OH' 
            %    'CH2CO'    'CH2CHO' 
            % 
            % indx = sol.lookForSpeciesPattern('CH2')
            % indx =
            %
            %    11    12    18    19    29    52
            %
            %             [index,species]  = b.lookForSpeciesPattern('C2','CH4','CH')
            %index =
            % 
            %   Columns 1 through 16
            % 
            %     10    11    12    13    14    18    19    20    21    22    23    24    25    26    27    29
            % 
            %   Columns 17 through 18
            % 
            %     52    53
            % 
            % 
            % species = 
            % 
            %   Columns 1 through 10
            % 
            %     'CH'    'CH2'    'CH2(S)'    'CH3'    'CH4'    'CH2O'    'CH2OH'    'CH3O'    'CH3OH'    'C2H'
            % 
            %   Columns 11 through 18
            % 
            %     'C2H2'    'C2H3'    'C2H4'    'C2H5'    'C2H6'    'CH2CO'    'CH2CHO'    'CH3CHO'
            % 
                        
            if nargin == 1
                obj.wrongNumberInputException.throwAsCaller;
            else
                pattern  = varargin(:);                       
            end  
            splistIn = obj.speciesNames;
            l = 0;
            indexList = [];
            if isa(pattern,'char')
                for kk = 1:length(splistIn)
                    if ~isempty(strfind(splistIn{kk},pattern))
                        l = l + 1;                        
                        indexList(l) = kk;
                    end
                end
            elseif isa(pattern,'cell')
                for k = 1:length(pattern)
                    for kk = 1:length(splistIn)                        
                        if ~isempty(strfind(splistIn{kk},pattern{k}))
                            l = l + 1;                            
                            indexList(l) = kk;
                        end
                    end
                end
            else
                error('findSpecies:WrongArgumentType',...
                    'The argument should be a string or a cell array of strings')
            end
            
            indexList =  unique(indexList);
            
            splistOut = obj.speciesNames(indexList);
            switch nargout
                case 0
                    
                    fprintf('\nFound ');
                    fprintf(num2str(length(indexList)));
                    fprintf(' species that contains ''');
                    fprintf(char(pattern));
                    fprintf('''\n');
                    chemistry.listTheSpecies(splistOut)
                case 1    
                    varargout{1} =  indexList;
                case 2
                    varargout{2} = splistOut;
                    varargout{1} = indexList;
                otherwise
            end       
        end        
         
        function varargout = mainSpecies(obj,varargin)
            % gives back the n main species of a system
            % main species is defined by the order of maximal mass/mole fraction
            % If a additional parameter n is given, n species are given back,
            % otherwise n = 6.
            % (C) 2011 U.Pruefert for VIRTUHCON
            
            % compute the main species
            % case 1: reactor not running...
            if isempty(obj.massFractions)
                if ~isempty(obj.initialMassFractions)
                    indx = obj.getSpeciesIndexByName(...
                        obj.getSpeciesFromInitialFractions(...
                        obj.initialMassFractions));
                      [~,indxRel] = sort(obj.getFractionsFromInitialFractions(...
                        obj.initialMassFractions),'descend');
                    indx = indx(indxRel(1:min(length(indx),6)));
                elseif ~isempty(obj.initialMoleFractions)
                    indx = obj.getSpeciesIndexByName(...
                        obj.getSpeciesFromInitialFractions(...
                        obj.initialMoleFractions));
                      [~,indxRel] = sort(obj.getFractionsFromInitialFractions(...
                        obj.initialMoleFractions),'descend');
                    indx = indx(indxRel(1:min(length(indx),6)));
                else
                    ME = MException('chemistry:EmptyData',...
                        'The object has empty fraction properties.');
                    throw(ME);
                end
            % case 2, mass/mole Fractions not empty
            else
                switch nargin
                    case 1
                        N = 6;
                        [~,indx] = sort(max(obj.moleFractions,[],2),'descend');
                        indx = indx(1:N);
                    case 2 % can be N or the mole/mass switch
                        switch class(varargin{1})
                            case 'char' 
                                N = 6;
                                switch varargin{1}                       
                                    case 'mole' 
                                        [~,indx] = sort(max(obj.moleFractions,[],2),'descend');
                                        indx = indx(1:N);                            
                                    case 'mass'
                                        [~,indx] = sort(max(obj.massFractions,[],2),'descend');
                                        indx = indx(1:N);    
                                    otherwise
                                        throwAsCaller(obj.wrongInputArgumentException);
                                end
                            case 'double'
                                N = varargin{1};
                                if N<1||N>obj.nSpecies||~isnumeric(N)
                                    error('chemistry:mainSpecies:checkInput',...
                                        'The number of main species must be between 1 and # species');
                                else
                                    [~,indx] = sort(max(obj.moleFractions,[],2),'descend');
                                    indx = indx(1:N);
                                end
                            otherwise
                        end
                    case 3
                        switch class(varargin{1})
                            case 'char'
                                type = varargin{1};
                            case 'double'
                                N = varargin{1};
                            otherwise
                        end
                        switch class(varargin{2})
                            case 'char'
                                type = varargin{2};
                            case 'double'
                                N = varargin{2};
                            otherwise
                        end
                        if N<1||N>obj.nSpecies||~isnumeric(N)
                            ME = MExcetption('chemistry:mainSpecies:checkInput',...
                                'The number of main species must be between 1 and # species');
                            throw(ME);
                        end
                        switch type
                            case 'mole' 
                                [~,indx] = sort(max(obj.moleFractions,[],2),'descend');
                                indx = indx(1:N);                            
                            case 'mass'
                                [~,indx] = sort(max(obj.massFractions,[],2),'descend');
                                indx = indx(1:N);    
                            otherwise
                                ME = MExcetption(...
                                    'chemistry:mainSpecies:checkInput',...
                                    'Second argument must be ''mole'' or ''mass''.');
                                throw(ME);
                        end
                    otherwise                        
                         throwAsCaller(obj.wrongNumberInputException);
                end
            end
            % return the result(s)
             switch nargout        
                 case 0
                     fprintf('The main species are: \n')
                     chemistry.listTheSpecies(obj.speciesNames(indx))
                 case 1
                     varargout{1} = indx';
                 case 2
                     varargout{1} = indx';
                     varargout{2} =  obj.speciesNames(indx);
                     
                 otherwise
             end
        end
        
        % convert mass to mole
        function moleFract = mole2mass(obj,moleFrac)
            % moleFractions = obj.mole2mass(Frac)
            % converts mole fractions to mass fractions in
            % cantera string format            
            % call some friends
            sp = obj.getSpeciesFromInitialFractions(moleFrac);
            mole = obj.getFractionsFromInitialFractions(moleFrac);
            % compute the index vector
            idx = obj.getSpeciesIndexByName(sp);
            mw = molecularWeights(obj.gas);
            mw = mw(idx);
            mass = mole.*mw' ;
            mass = mass./sum(mass);
            moleFract  =  obj.vector2CanteraString(sp,mass);
        end
        
        function moleFract = mass2mole(obj,massFrac)
            % moleFractions = obj.mass2mole(massFrac)
            % converts massfractions to mole fractions in
            % cantera string format
            
            % call some friends
            sp = obj.getSpeciesFromInitialFractions(massFrac);
            mass = obj.getFractionsFromInitialFractions(massFrac);
            % compute the index vector
            idx = obj.getSpeciesIndexByName(sp);
            mw = molecularWeights(obj.gas);
            mw = mw(idx);
            moles = mass./mw' ;
            moles = moles./sum(moles);
            moleFract  =  obj.vector2CanteraString(sp,moles);
        end        
         
        function   removeSpecies(obj,varargin)
            %  removeSpecies(mechFileOut,deadSpeciesIndexList)
            %  removeSpecies(baseMech,mechFileOut,deadSpeciesIndexList,deadSpeciesNamesList)
            %   removes species given by deadSpeciesIndexList and deadSpeciesNamesList
            %   from a mechanism file. If baseMech is a cti file, removeSpeciesCTI
            %   will used, is mechFIleIn a inp-file, removeSpeciesCK will be used. The
            %   type of the format of the input file is detected by the ending inp or
            %   cti.
            %   New Version.
            %   (C) 2013 Uwe Prüfert 
            
            switch nargin
                case 3
                    baseMech = obj.mechanism;        
                    mechFileOut = varargin{1};
                    if iscell(varargin{2})
                        deadSpeciesNamesList = varargin{2};
                        deadSpeciesIndexList = obj.getSpeciesIndexByName(...
                            deadSpeciesNamesList);
                    else
                        deadSpeciesIndexList = varargin{2};
                        deadSpeciesNamesList = obj.speciesNames(...
                            deadSpeciesIndexList);
                    end
                case 4
                    baseMech =  varargin{1};       
                    mechFileOut = varargin{2};
                    deadSpeciesNamesList = varargin{3};
                    if ~strcmp(baseMech(end-2:end),'inp')
                        fprintf('three arguments call is only possible for Chemkin-files.')
                        obj.wrongInputClassException.throw
                    end
                case 5
                    baseMech =  varargin{1};       
                    mechFileOut = varargin{2};
                    deadSpeciesIndexList = varargin{3};
                    deadSpeciesNamesList = varargin{4};
                otherwise
            end
            % identify the input file type
            dotloc = max(strfind(mechFileOut,'.'));
            if dotloc > 1
                filetype = mechFileOut(dotloc+1:end);
            end
            
            switch filetype
                
                case 'inp'
                    fprintf('Create Chemkin file\n');
                    % call CK version
                    obj.removeSpeciesCK(baseMech,mechFileOut,...
                        deadSpeciesNamesList)
                case 'cti'
                    fprintf('Create Cantera file\n');
                    % call CTI version  
                    obj.removeSpeciesCTI(baseMech,mechFileOut,...
                        deadSpeciesIndexList,deadSpeciesNamesList)
                otherwise
                    error('removeSpecies:wrongFileType',...
                        'Input-Output can be *.inp or *.cti')
            end            
        end
        
        
    end   % public methods
      
    methods(Static,Access = public)
        function [varargout] = getSpeciesFromMech(mechfile)
            %[speciesList[mechFile]]] = getSpeciesFromMech(mechfile)
            %   reads the species from a *.cti or *.inp file. The optional output argument contains the
            %   the path to the mechfile.
            nargs = nargout;  
            if nargin == 0
                throwAsCaller(chemistry.wrongNumberInputException);
            end
            fid = fopen(mechfile,'r');
            if fid<0
                 throwAsCaller(chemistry.noFileException)
            end
            dotloc = max(strfind(mechfile,'.'));
            if dotloc > 1
                filetype = mechfile(dotloc+1:end);
            else
                throwAsCaller(chemistry.wrongInputFileFormatException)                
            end
            switch filetype
                case 'inp'
                    k = 0;
                    nspecies = 0;
                    speciesList = []; 
                    inThermo = false;
                    % call CK version
                    while 1
                        cline = fgetl(fid); 
                        if ~ischar(cline),         
                            break,
                        end
                        if (length(cline)>=3) && (strcmp(cline(1:3),'END')||strcmp(cline(1:3),'end'))
                            inThermo = false;   % exit all blocks
                        end
                        if inThermo
                            % in thermo blocks when a line starting with a letter, the
                            % fist charakters are a definition of a species name
                            id = double(cline(1));
                            if (id>=65&&id<=90)||(id>=97&&id<=122)
                                % it must be a letter hence it is a NAME
                                endStr = strfind(cline,' ');
                                k = k+1;
                                nspecies = nspecies + 1;
                                % the fist blank
                                speciesList{k} = cline(1:endStr-1);
                            end
                        end
                        if (length(cline)>=6) && (strcmp(cline(1:6),'THERMO')||strcmp(cline(1:6),'thermo'))
                            % reactions block
                            inThermo = true;
                        end
                    end
                case 'cti' % handle CTI file
                    k = 1;
                    nspecies = 0;
                    speciesList = [];        
                    while 1
                        tline = fgetl(fid);
                        if length(tline)>16 && strcmp(tline(1:7),'species')
                            nspecies =  nspecies + 1;
                            species = [];
                            for l = 17:80,
                                if strcmp(tline(l),'"')
                                    break
                                else
                                    species = [species,tline(l)];
                                end
                            end
                            speciesList{k} = species;
                            k = k+1;    
                        end
                        if ~ischar(tline)
                            % return statement
                            break
                        end
                    end
                otherwise
                    throwAsCaller(obj.wrongInputFileFormatException);
            end
            if isempty(speciesList)
                warning('getSpeciesFromMech:EMptyList',...
                    'The species list is empty, Maybe you have a corrupted or empty mechanism file.')
            end
            switch nargs
                case 1
                    varargout{1} = speciesList;
                case 2
                    varargout{1} = speciesList;
                    varargout{2} = mechfile;
                otherwise
            end
        end   
             
        function varargout = isInList(list,str)
            % bool = isInList(list,str) returns true if the 
            % string str ist found in the cell array list. Case sentitive! 
            b = false;
            for kk = 1:length(list)
                b = strcmp(list{kk},str);
                if b        
                    break
                end
            end
            switch nargout
                case 1
                    varargout{1} = b;
                case 2
                    varargout{1} = b;
                    varargout{2} = kk;
                otherwise
                    varargout = [];
            end
        end
        
        function [b,varargout] = isInListI(list,str)
            % bool = isInListI(list,str) returns true if the 
            % string str ist found in the cell array list. NOT case sensitive!!!
            b = false;
            for kk = 1:length(list)
                b = strcmpi(list{kk},str);
                if b        
                    break
                end
            end
            switch nargout
                case {0 1}
                    % b = b :-)
                case 2       
                    varargout{1} = kk;
                otherwise
                    MException('chemistry:isInListI:WrongNumberOutputs',...
                        'The number of outputs can be one or two.').throw;
            end
        end
       
        function fractionsOut= mixFractions(fractions1,fractions2,a)
            % This function mix the fractions given by the firts two arguments,
            % asssumed as a valid fractions strings.
            % The third argument is the  ratio phi. We put 1/phi of the mixture in
            % the first argument to the initial mass fractions in
            % obj.initialMassFraction
            % (C) 2011 by Uwe Pruefert for VIRTUHCON
            if a>1 || a <0
                error('matlab:batchReactor:friends:mixFractions',...
                    'Wrong value for mixing ratio: must be between 0 and 1.')
            end
            speciesListOne = chemistry.getSpeciesFromInitialFractions(fractions1);
            fractionsListOne = a * chemistry.getMoleFractionsFromInitialFractions(fractions1)/...
                sum(chemistry.getMoleFractionsFromInitialFractions(fractions1));
            
            speciesListTwo = chemistry.getSpeciesFromInitialFractions(fractions2);
            fractionsListTwo =  (1-a) * chemistry.getMoleFractionsFromInitialFractions(fractions2)/...
                sum(chemistry.getMoleFractionsFromInitialFractions(fractions2));
            
            % mix the species and the fractions
            
            newSpecies = logical(1:length(speciesListTwo));
            
            for k = 1:length(speciesListOne)    
                for l = 1:length(speciesListTwo)        
                    if strcmp(speciesListOne{k},speciesListTwo{l})            
                        fractionsListOne(k) =  fractionsListOne(k) + fractionsListTwo(l);
                        newSpecies(l) = false;           
                        break;        
                    end
                end
            end
            
            % new species is in the end a logical vector where true is when the species
            % name is not fpount in the initial mass fractions of batchReactor object
            
            speciesListOne = [speciesListOne speciesListTwo(newSpecies)];
            fractionsListOne = [fractionsListOne fractionsListTwo(newSpecies)];
            
            fractionsOut  = chemistry.vector2CanteraString(...
                speciesListOne,fractionsListOne);
        end
            
        function blocks = indexCtiFile(ctiFile)
            % blocks = indexCtiFile(ctiFile)
            % indexes a cti file
            % 
            % blocks.nlines  : no of lines;
            % blocks.species : startline of species block
            % blocks.ideal_gas :  startline ideal_gas_block
            % blocks.reactions :  startline of a reaction_block
            % blocks.three_bodies :  startline of thethree_body_reaction_block
            % blocks.falloffs     :  startlinefalloff_reaction_block
            % blocks.units        :  startline of the unit_block            
            fid = fopen(ctiFile,'r');
            nline = 0;
            ideal_gas_block = [];
            species_block = [];
            reaction_block = [];
            falloff_reaction_block = [];
            three_body_reaction_block = [];
            unit_block = [];
            reaction_index = [];
            while 1
                cline = fgetl(fid);     %read line
                nline = nline+1;        % count line                
                % identify the keywords and the associate blocks
                if length(cline)>=10 && strcmp(cline(1:9),'ideal_gas')
                    ideal_gas_block = [ideal_gas_block,nline];
                elseif length(cline)>=7 && strcmp(cline(1:7),'species')
                    species_block = [species_block,nline];
                elseif length(cline)>=8 && strcmp(cline(1:8),'reaction')
                    reaction_block = [reaction_block,nline];
                    reaction_index = [reaction_index,nline];
                elseif length(cline)>=16 && strcmp(cline(1:16),'falloff_reaction')
                    falloff_reaction_block = [falloff_reaction_block,nline]; %#ok<*AGROW>
                    reaction_index = [reaction_index,nline];
                elseif length(cline)>=19 && strcmp(cline(1:19),'three_body_reaction')
                    three_body_reaction_block = [three_body_reaction_block,nline];
                    reaction_index = [reaction_index,nline];
                elseif length(cline)>=5 && strcmp(cline(1:5),'units')
                    unit_block = nline;
                end                
                if ~ischar(cline),
                    err = fclose(fid);
                    break,
                end
            end
            blocks.nlines = nline-1;
            blocks.species = species_block;
            blocks.nspecies = length(blocks.species);
            blocks.ideal_gas = ideal_gas_block;
            blocks.reaction_index = reaction_index; 
            blocks.reactions = reaction_block;
            blocks.three_bodies = three_body_reaction_block;
            blocks.falloffs = falloff_reaction_block;
            blocks.units = unit_block;
            blocks.nreactions = length(blocks.reactions)...
                +length(blocks.three_bodies)+length(blocks.falloffs);            
        end
        
        function CanteraString = vector2CanteraString(speciesNames,v)
            % CanteraString = vector2CanteraString(speciesNames,vectorOfFractions)
            % Function to create a Mass/Molefractions string for use with canteras
            % set(gas) method. Useful to initialize a gas object with values computed by
            % another program or given data.
            % Exapmle:
            % Assume we have a chemistry object c with a solved reactor
            % We want to use send the mole fraction at the end of the time domain to a
            % gas object.
            %   moles = c.getMoleFractions({'CH4' 'O2'});
            %   set(gas,'X',vector2CanteraString({'CH4' 'O2'},moles(:,end)));
            %  (c) 2013 by   Uwe Prüfert 
            if ~iscellstr(speciesNames)
                error('Vector2CanteraString:MustBeCell','The Names of the species must be a cell array of strings');
            end
            if ~isnumeric(v)
                error('Vector2CanteraString:NotNumeric','The array of the fractions must be numeric');
            else
                if length(speciesNames)~=length(v)
                    error('Vector2CanteraString:DimensionMissmatch','The dimensions of the names and the fractions must be equal');
                end
            end
            CanteraString = '';
            n = length(speciesNames);
            for i = 1:n
                S = char(speciesNames(i));
                CanteraString = strcat(CanteraString,S,':',num2str(v(i)));
                if i~=n
                    CanteraString=strcat(CanteraString,',');
                end
            end
        end
        
    end
        
    methods(Static,Access = protected)
        
        function newEfficiencies = removeDeadFromEfficiencies(oldEfficiencies,deads)
            % scans efficiencies list and remove species in dead
            l = 1;m=1;remove = false;
            startName = false;
            for k = 26:length(oldEfficiencies)
                if isLetter(oldEfficiencies(k)) && ~startName
                    startSpeciesName(l) = k;
                    startName = true;
                    l = l+1;
                elseif strcmp(oldEfficiencies(k),':')
                    endSpeciesName(m) = k-1;
                    startName = false;
                    m = m+1;
                end
            end
            endSpecies = startSpeciesName(2:end)-3;
            endSpecies(length(startSpeciesName))=length(oldEfficiencies)-3;

            if ~(length(startSpeciesName)==length(endSpeciesName))
                MException('chemistry:WrongLengthArguments',...
                    'The format of the input file may be corrupted.').throw;
            end
            newEfficiencies = oldEfficiencies(1:25) ;
            for k = 1:length(startSpeciesName),
                for l = 1:length(deads)  
                    if strcmp(oldEfficiencies(startSpeciesName(k):endSpeciesName(k)),deads{l}),
                        remove = true;           
                        break
                    else
                        remove = false;
                    end
                end
                if ~remove
                    newEfficiencies = [newEfficiencies,' ',oldEfficiencies(startSpeciesName(k):endSpecies(k)),' ' ];
                end
            end
            newEfficiencies = [newEfficiencies,'")']  ;
            function b = isLetter(charSpec)
                % is it a  letter character (to define a species name, the first must
                % be a letter, 20H is NOT a well formed species name
                b = (double(charSpec)<91 && double(charSpec)>64) ||...% capitals A...Z
                            (double(charSpec)<123 && double(charSpec)>96);  % lowercasses a..z
            end
        end

        function speciesList = getSpeciesFromInitialFractions(str) 
            %speciesList = getSpeciesFromInitialFractions(str) 
            % Extracts species names from the (user input) inital mole fractions
            % here we use it for selecting the three body from their definiton
            % sequence what is quite similar to the cantera systax of the mole
            % fractions.
            % Note, that "," is NOT allowed. However, it can be a part of some
            % strange species names...  
             
            if isempty(str)
                throwAsCaller(obj.noStringException);
            end
            speciesList = [];
            k = 1; % zaehler fuer die position
            in = true;
            species = [];
            for i=1:length(str)        
                if strcmp(str(i),':')
                    in = false; 
                    speciesList{k} = species;
                    species = [];
                    k = k+1; % next  
                end
                if strcmp(str(max(1,i-1)),',')||strcmp(str(max(1,i-1)),'')
                    in = true;
                end
                if allowedChar(str(i)) && in
                    species = [species,str(i)];             
                end
            end
            % local ++++
            function b = allowedChar(charSpec)
                % is it a allowed character in the definietion of a species?   
                b = (double(charSpec)<91 && double(charSpec)>64) ||...% capitals A...Z
                    (double(charSpec)<123 && double(charSpec)>96) ||...lowercasses a..z
                    (double(charSpec)<59 && double(charSpec)>47) ||... % numbers 0,1...9
                    ...  %double(charSpec)==44||... % ,
                    double(charSpec)==95||... % _
                    double(charSpec)==35||... % -
                    double(charSpec)==45||... % -
                    double(charSpec)==43||... % +
                    double(charSpec)==41||... % (
                    double(charSpec)==40;     % )                 
            end
        end
                
        function moles = getFractionsFromInitialFractions(str,varargin)
            % extracts molefractions from the (user input) inital mole
            % fractions and normalze them (if necessary)
            if isempty(str)
                throwAsCaller(obj.noStringException);
            end
            speciesMole =[];
            k = 1; 
            in = false;
            for i=1:length(str)   
                if strcmp(str(i),':')
                    in = true;               
                end
                if isnumber(str(i)) && in
                    speciesMole = [speciesMole,str(i)];             
                end
                if (strcmp(str(i),',')||i==length(str)||strcmp(str(i),' ')) && in
                    moles(k) = str2double(speciesMole);
                    k=k+1;
                    in = false;
                    speciesMole=[];
                end
            end
            switch nargin
                case 3
                    if strcmp(varargin{1},'normalize')
                        if strcmp(varargin{2},'off')  
                            %                
                        elseif strcmp(varargin{2},'on')
                            moles = moles./sum(moles);
                        end
                    end
                otherwise
                    moles = moles./sum(moles);
            end
            % local
            function b = isnumber(char)
                % checks if char is a number, i.e. it contains '0'...'9' and '.'
                b = ((double(char)<58 && double(char)>47))||double(char)==46||double(char)==45||double(char)==101;
            end
        end 
        
        function listTheSpecies(list)
            m = ceil(length(list)/4);
            r = length(list);
            n = 1;
            fprintf('\n')
            for k = 1:m,
                for l = 1:4,
                    if n <= r
                        fprintf(['\t',list{n}]);
                        n = n+1;
                    end
                end
                fprintf('\n')
            end
            fprintf('\n')
        end
             
        function nline = removeDeadSpeciesFromIdealGasBlock(reactionLine,deads)
            %nline = removeDeadSpeciesFromIdealGasBlock(reactionLine,deads)
            %   We remove the species from the dead list and write back a new species
            %   definition line

            % extract the species from the line after the keyword
            species = chemistry.scanSpeciesList(reactionLine(20:end));

            newSpecies=[];


            for l = 1:length(species)
                for m = 1:length(deads)
                    if strcmp(species{l},deads{m})
                       species{l} = '';         
                    else
                        %fprintf(['remove species ',species{l},' from list\n'])
                    end
                end
            end
            l = 1;
            for k = 1:length(species)
                if ~isempty(species{k})
                    newSpecies{l}=species{k};
                    l = l+1;
                end
            end
            newSpeciesString = [];
            for k = 1:length(newSpecies)
                newSpeciesString = [newSpeciesString,' ',newSpecies{k},' '];
            end
            % we add the new line to the first 19 entries in the old line, i.e we write
            % " species" keyword and preserve the structure of the file...
            % we also have to set the block closing '""",'
            if strcmp(reactionLine(end-3:end),'""",')
                nline = [reactionLine(1:18),newSpeciesString,'""",'];
            elseif strcmp(reactionLine(7:13),'species')

                nline = [reactionLine(1:19),newSpeciesString];
            else

                nline = [reactionLine(1:19),newSpeciesString];
            end
        end
        
        function removeReactionsCK(mechFileIn,mechFileOut,varargin)
        % removeReactionsCK(mechFileIn,mechFileOut[,species_list])
        % remove all reactions that contains the species in species_list.
        % species_list must be a CELL array, e.g. {'H2' 'CO2'} or a species object.
        % If no species_list is given, removeReactions removes ALL reactions from 
        % the mechanism.
        % (c) U. Pruefert (VIRTUHCON) 2011.

        % check for stupid inputs/output file declaration
        if strcmp(mechFileIn,mechFileOut)
            error('removeDeadSpecies:SameInputOutputFile',...
                'The input and the output files are the same, are you sure? However, it is not allowed. Bye!')
        end

        % try to open the files
        [fidMech,mesg] = fopen(mechFileIn,'r');
        if fidMech < 0 
            fprintf(mesg)
        end 
        [fidNewMech,mesg] = fopen(mechFileOut,'w');
        if fidNewMech < 0 
            fprintf(mesg)
        end

        % some comment add to the mechanism file
        fprintf(fidNewMech,['! This is a reduced mechanism based on ',mechFileIn,' \n']);
        fprintf(fidNewMech,'! The removing process is performed by removeReactionsCK\n');
        fprintf(fidNewMech,'! \t\t (c) U. Pruefert 2011.\n!\n');

        if isempty(varargin)
            all = true;
        else
            if isa(varargin{1},'cell')
                deadSpeciesNamesList = varargin{1}; 
            elseif isa(varargin{1},'species')
                deadSpeciesNamesList = varargin{1}.names;
            else
                error('removeReactionsCK:wrongArgument',...
                    'The optinal argument must be a CELL array containing the names of species or a SPECIES object')
            end
            all = false;
            fprintf(fidNewMech,'! All reactions associated with the following species will be deleted:\n');
            listTheSpecies(fidNewMech,deadSpeciesNamesList)
        end
        inReactions = false;
        toRemove = true;
        nOfReactions = 0;nOfRemovedReactions = 0;
        while 1 % loop over all lines
            cline = fgetl(fidMech); 
            if ~ischar(cline),         
                break,
            end 
            if (length(cline)>=4) && (strcmp(cline(1:4),'REAC')||strcmp(cline(1:4),'reac'))
                % reactions block
                inReactions = true;
                fprintf(fidNewMech,[cline,'\n']);
            elseif (length(cline)>=3) && (strcmp(cline(1:3),'END')||strcmp(cline(1:3),'end')) 
                % exit all blocks         
                inReactions = false;

                fprintf(fidNewMech,[cline,'\n']);
            else
                % now, we have to work...
                % Definitely, we are in a non keyword line or a (regular) comment
                % line
                if inReactions
                    % first we look for a line that defines a reaction.
                    % if a reaction contains a to-removed-species, all lines can be
                    % removed untill we found the next reaction
                    if isReactionLine(cline)
                        if ~all
                            toRemove = false;
                        end
                        nOfReactions = nOfReactions + 1;
                    end
                    if ~toRemove
                        if containsDeadSpecies(cline,deadSpeciesNamesList)
                            toRemove = true;
                            nOfRemovedReactions = nOfRemovedReactions + 1;
                        end
                    end
                    if ~toRemove                
                        fprintf(fidNewMech,[cline,'\n']);
                    end
                else
                    fprintf(fidNewMech,[cline,'\n']);            
                end
            end   
        end

        % we close all files
        err = fclose(fidNewMech);
        if err < 0
            fprintf('An error during closing the file is occured. Please check the file names etc.\n')
        end
        err = fclose(fidMech);
        if err < 0
            fprintf('An error during closing the file is occured. Please check the file names etc.\n')
        end



            function b = isReactionLine(cline)
                b = false;
                for l = 1:length(cline)-1;
                    % if we found a comment mark, we can abort the procedure
                    if strcmp(cline(l),'!')
                        return
                    end
                    % a reaction line contains "=" (within "=>", "<=", "<=>", or "=" )
                    b = strcmp(cline(l),'=');
                    if b 
                        return                          
                    end
                end
            end

        % we check, if a line contains a species
            function b = containsDeadSpecies(cline,deadSpeciesNamesList)
                b = false;
                for k = 1:length(deadSpeciesNamesList)
                    lngthSpecName = length(deadSpeciesNamesList{k});
                    for l = 1:length(cline)-lngthSpecName-1;
                        foundSpec = strcmp(cline(l:l+lngthSpecName-1),deadSpeciesNamesList{k});
                        b = isLeft(cline,l) && isRight(cline,l) && foundSpec ; 
                        if b 
                            return
                        end               
                    end
                end 
                % local local functions for identifying species by name
                function b = isLeft(cline,l)
                    % checks if a species string is a whole species name or only a part
                    % of a species name: left of the name...
                    b = (l==1) || ((l>1) && (strcmp(cline(l-1),'+') ||  strcmp(cline(l-1),'>')||...
                        strcmp(cline(l-1),'=') || (double(cline(l-1))>48 && double(cline(l-1))< 58) )...
                         &&...
                        ((l>=2)||(strcmp(cline(l-2),'+') ||  strcmp(cline(l-2),'>') || strcmp(cline(l-2),'=')))); 
                end
                % of a species name: right of the name...
                function b = isRight(cline,l)
                    b = ( strcmp(cline(l+lngthSpecName),' ') ||  strcmp(cline(l+lngthSpecName),'<') ||...
                        strcmp(cline(l+lngthSpecName),'+')||  strcmp(cline(l+lngthSpecName),'=') ||...
                        strcmp(cline(l+lngthSpecName:l+lngthSpecName+1),'(+')) ; 
                end
            end

        % helper function that lists the species
            function listTheSpecies(fid,list)
                m = ceil(length(list)/4);
                r = length(list);
                n = 1;
                fprintf(fid,'!\n');
                for k = 1:m
                    fprintf(fid,'!');
                    for l = 1:4,
                        if n <= r
                            fprintf(fid,['\t',list{n}]);
                            n = n+1;
                        end
                    end
                    fprintf(fid,'\n');
                end
                fprintf(fid,'! \n');
            end
        end
        
        function removeReactionsCTI(mechFileIn,mechFileOut,varargin)

            if strcmp(mechFileIn,mechFileOut)
                error('removeReactionsCTI:SameInputOutputFile','The input and the output files are the same, are you sure? However, it is not allowed. Bye!')
            end

            blocks = indexCtiFile(mechFileIn);

            % Check whether any index array is empty
            % If there are no reactions, three-bodies or falloffs, we set the line to 1e5.
            % Then the mechanism is already a non-reaction mechanism and we will copy
            % it. reactionsStartAt is now larger then the number of lines.
            if isempty(blocks.reactions)
                reactions = 1e5;
            else
                reactions = blocks.reactions;
            end
            if isempty(blocks.three_bodies)
                three_bodies = 1e5;
            else
                three_bodies = blocks.three_bodies;
            end
            if isempty(blocks.falloffs)
                falloffs = 1e5;    
            else
                falloffs = blocks.falloffs;
            end

            % a boolean for swithcind all <=> single reaction remove
            removeSingleReactions = nargin>2;
            isToRemove = [];
            if removeSingleReactions
                indexSet = varargin{1};
                % create the vector of lines to be removed
                allReactions = sort([blocks.reactions blocks.three_bodies blocks.falloffs]);

                for k = 1:length(indexSet)
                    if indexSet(k)+1>length(allReactions)
                        isToRemove = [isToRemove allReactions(indexSet(k)):allReactions(indexSet(k))+6];
                    else
                        isToRemove = [isToRemove allReactions(indexSet(k)):allReactions(indexSet(k)+1)-1];
                    end
                end
            end



            reactionsStartAt = min([reactions(1) ...
                three_bodies(1) falloffs(1) ]);

            % tries to open the files
            [fidMech,mesg] = fopen(mechFileIn,'r');
            if fidMech < 0 
                fprintf(mesg)
            end 
            [fidNewMech,mesg] = fopen(mechFileOut,'w');
            if fidNewMech < 0 
                fprintf(mesg)
            end

            % some comment add to the mechanism file
            fprintf(fidNewMech,['# This is a reduced mechanism based on ',mechFileIn,' \n']);
            fprintf(fidNewMech,'# The removing process is performed by removeReactionsCTI\n');
            fprintf(fidNewMech,'# \t\t (c) U. Pruefert 2011.\n\n');

            nLine = 0;

            while 1 % loop over all lines
                cline = fgetl(fidMech); 
                nLine = nLine + 1;
                if ~ischar(cline)
                    break
                end
                if (~ischar(cline)||nLine==reactionsStartAt-1)&&~removeSingleReactions,         
                    break,
                elseif isIn(nLine,isToRemove)
                    % do nothing...
                else
                    fprintf(fidNewMech,[cline,'\n']);
                end
            end

            % we close all files
            err = fclose(fidNewMech);
            if err < 0
                fprintf('An error during closing the file is occured. Please check the file names etc.\n')
            end
            err = fclose(fidMech);
            if err < 0
                fprintf('An error during closing the file is occured. Please check the file names etc.\n')
            end

                function b = isIn(x,y)
                    % true if x is in y
                    if isempty(y)
                        b = false;
                        return
                    end
                    for kk = 1:length(y)
                        b = x==y(kk);
                        if b
                            break
                        end
                    end
                end

            end
             
        function b = isDeadReaction(reactionLine,deads)
            % b= isDeadReaction(reactionLine) 
            %   decides if an reaction is dead, i.e. if it contains a dead species
            b = false;
            species = chemistry.selectSpeciesFromList(reactionLine);
            for k = 1:length(deads)
                for l = 1:length(species)
                    if strcmp(deads{k},species{l})
                        b = true;
                        return
                    end
                end
            end
        end
     
        function removeSpeciesCK(mechFileIn,mechFileOut,varargin)
            %removeSpeciesCK(mechFileIn,mechFileOut,deadSpeciesNamesList)
            % This is the CHEMKIN Version of removeSpecies.
            %   Removes the species (and Reactions) defined in 
            %   deadSpeciesNamesList from mechFileIn, the resulting mechanism will be
            %   written in mechFileOut. 
            %   The input arguments deadSpeciesNamesList must be
            %   either a CELL array containing the names of the species, e.g. {'OH' 'CH3'}
            %   or a species object  
            %   (c) U. Pruefert (VIRTUHCON) 2011. 
            % check for stupid inputs/output file declaration
            if strcmp(mechFileIn,mechFileOut)
                MException('removeDeadSpecies:SameInputOutputFile',...
                    ['The input and the output files are the same,',...
                    ' are you sure? However, it is not allowed. Bye!']).throw;
            end
            
            if isa(varargin{1},'cell')
                deadSpeciesNamesList = varargin{1}; 
            elseif isa(varargin{1},'species')
                deadSpeciesNamesList = varargin{1}.names;
            else
                MException('removeSpeciesCK:wrongArgument',...
                    ['The optinal argument must be a CELL',...
                    ' array containing the names of species or a SPECIES object']).throw;
            end
            % try to open the files
            [fidMech,mesg] = fopen(mechFileIn,'r');
            if fidMech < 0 
                fprintf(mesg)
            end
            [fidNewMech,mesg] = fopen(mechFileOut,'w');
            if fidNewMech < 0 
                fprintf(mesg)
            end
            
            % initalize the line counter
            nLines = 0;
            nOfReactions = 0;
            nOfRemovedReactions = 0;
            
            % some important counters
            inElements = false;
            inSpecies = false;
            inThermo = false;
            inReactions = false;
            linesJetToRemove = 0;
            
            % write some comments in the mechanism file
            fprintf(fidNewMech,['! This is a reduced mechanism based on ',mechFileIn,' \n']);
            fprintf(fidNewMech,'! The following species are deleted:\n');
            listTheSpecies(fidNewMech,deadSpeciesNamesList)
            fprintf(fidNewMech,'! Also, all reactions that contain species listed above are removed.\n');
            fprintf(fidNewMech,'! The removing process is performed by removeSpeciesCK\n');
            fprintf(fidNewMech,'! \t\t (c) U. Pruefert 2011.\n!\n');
            

            while 1 % loop over all lines
                cline = fgetl(fidMech); 
                % because we use a fixed condition in our while loop, we need a loop
                % break statement: if the next "line"contains NO character, it can only
                % be a EOF
                if ~ischar(cline),         
                    break,
                end
                nLines = nLines + 1;
                if (length(cline)>=1) && strcmp(cline(1),'!')
                    % comment line, write it in the output or not???
                    fprintf(fidNewMech,[cline,'\n']);
                elseif (length(cline)>=4) && (strcmp(cline(1:4),'ELEM')||strcmp(cline(1:4),'elem')) 
                    inElements = true;
                    fprintf(fidNewMech,[cline,'\n']);
                elseif (length(cline)>=4) && (strcmp(cline(1:4),'SPEC')||strcmp(cline(1:4),'spec')) 
                    % species block: Remove all species from the to-remove-list
                    inSpecies = true;
                    fprintf(fidNewMech,[cline,'\n']);
                elseif (length(cline)>=4) && (strcmp(cline(1:4),'THER')||strcmp(cline(1:4),'ther')) 
                    inThermo = true;
                    fprintf(fidNewMech,[cline,'\n']);
                elseif (length(cline)>=4) && (strcmp(cline(1:4),'REAC')||strcmp(cline(1:4),'reac'))
                    % reactions block
                    inReactions = true;
                    fprintf(fidNewMech,[cline,'\n']);
                elseif (length(cline)>=3) && (strcmp(cline(1:3),'END')||strcmp(cline(1:3),'end')) 
                    % exit al block
                    inSpecies = false;
                    inThermo = false;
                    inReactions = false;
                    inElements = false;
                    fprintf(fidNewMech,[cline,'\n']);
                else
                    % now, we have to work...
                    % Definitely, we are in a non keyword line or a (regular) comment
                    % line
                    if inElements 
                        % only to copy the element block 
                        fprintf(fidNewMech,[cline,'\n']);            
                    elseif inSpecies
                        % remove all species from the line and fill it with blanks
                        fprintf(fidNewMech,[removeSpeciesFromLine(cline,deadSpeciesNamesList),'\n']);   
                    elseif inThermo
                        % a thermo block has exactly 4 lines. If we identify a to-rmove
                        %-species in a block, the next three lines are also be removed
                        if thermoBlockToRemove(cline,deadSpeciesNamesList)
                            linesJetToRemove = 3;
                        elseif linesJetToRemove > 0               
                            linesJetToRemove = linesJetToRemove - 1;
                        else
                            fprintf(fidNewMech,[cline,'\n']);
                        end
                    elseif inReactions
                        % first we look for a line that defines a reaction.
                        % if a reaction contains a to-removed-species, all lines can be
                        % removed untill we found the next reaction
                        if isReactionLine(cline)
                            toRemove = false;
                            nOfReactions = nOfReactions + 1;
                        
                        end
                        if ~toRemove
                            if containsDeadSpecies(cline,deadSpeciesNamesList)
                                toRemove = true;
                                nOfRemovedReactions = nOfRemovedReactions + 1;
                            end
                        end
                        if ~toRemove
                            % maybe we must remove species from the efficiencies line
                            if isEfficienciesLine(cline) && efficienciesContainsDeadSpecies(cline,deadSpeciesNamesList)
                                cline = removeSpeciesFromEfficiencies(cline,deadSpeciesNamesList);                   
                            end
                            fprintf(fidNewMech,[cline,'\n']);
                        end
                    end
                end
            end
            
            % we close all files
            fclose('all');
             

            % some statitics
            fprintf('\tremoveSpeciesCK: I removed the following ') 
            fprintf([num2str(length(deadSpeciesNamesList)),' species:\n']);
            listTheSpecies(1,deadSpeciesNamesList)
            fprintf(['\n and\n ',num2str(nOfRemovedReactions),' of ',num2str(nOfReactions),' reactions.\n']);


            % And now: some internal functions
            % we replace the to-removed-species in the species block by blanks
                function cline = removeSpeciesFromLine(cline,deadSpeciesNamesList)
                    for k = 1:length(deadSpeciesNamesList)
                        lngthSpecName = length(deadSpeciesNamesList{k});
                        for l = 1:length(cline)-lngthSpecName + 1;                
                            foundSpec = strcmp(cline(l:l+lngthSpecName-1),deadSpeciesNamesList{k});
                            if (((l+lngthSpecName)>length(cline) ||  strcmp(cline(l+lngthSpecName),' ')) && foundSpec)...
                                    && ((l==1 ||  strcmp(cline(l-1),' ')) && foundSpec)
                                cline(l:l+lngthSpecName-1) = blank(lngthSpecName);
                            end                 
                        end
                    end        
                end

            % we identify the species in the thermo definition line
                function b = thermoBlockToRemove(cline,deadSpeciesNamesList)
                    b = false;
                    for k = 1:length(deadSpeciesNamesList)
                        lngthSpecName = length(deadSpeciesNamesList{k});
                        b = strcmp(cline(1:lngthSpecName),deadSpeciesNamesList{k})&&strcmp(cline(lngthSpecName+1),' ');
                        if b               
                            return
                        end

                    end         
                end

            % we check if a line is a reaction definition line
                function b = isReactionLine(cline)
                    b = false;
                    for l = 1:length(cline)-1;
                        % if we found a comment mark, we can abort the procedure
                        if strcmp(cline(l),'!')
                            return
                        end
                        % a reaction line contains "=" (within "=>", "<=", "<=>", or "=" )
                        b = strcmp(cline(l),'=');
                        if b 
                            return                          
                        end
                    end
                end

            % we check, if a line contains a species
                function b = containsDeadSpecies(cline,deadSpeciesNamesList)
                    b = false;
                    for k = 1:length(deadSpeciesNamesList)
                        lngthSpecName = length(deadSpeciesNamesList{k});
                        for l = 1:length(cline)-lngthSpecName-1;
                            foundSpec = strcmp(cline(l:l+lngthSpecName-1),deadSpeciesNamesList{k});
                            b = isLeft(cline,l) && isRight(cline,l) && foundSpec ; 
                            if b 
                                return
                            end               
                        end
                    end 
                    % local functions for identifying species by name
                    function b = isLeft(cline,l)
                        % checks if a species string is a whole species name or only a part
                        % of a species name: left of the name...
                        b = (l==1) || ((l>1) && (strcmp(cline(l-1),'+') ||  strcmp(cline(l-1),'>')||...
                            strcmp(cline(l-1),'=') || (double(cline(l-1))>48 && double(cline(l-1))< 58) )...
                             &&...
                            ((l>=2)||(strcmp(cline(l-2),'+') ||  strcmp(cline(l-2),'>') || strcmp(cline(l-2),'=')))); 
                    end
                    % of a species name: right of the name can be a blank, "<", "+" "=" "(+", or "/"...
                    function b = isRight(cline,l)
                        b = ( strcmp(cline(l+lngthSpecName),' ') ||  strcmp(cline(l+lngthSpecName),'<') ||...
                            strcmp(cline(l+lngthSpecName),'+')||  strcmp(cline(l+lngthSpecName),'=') ||...
                            strcmp(cline(l+lngthSpecName:l+lngthSpecName+1),'(+') ||...
                            strcmp(cline(l+lngthSpecName),'/')) ; 
                    end
                end

        % helper function that lists the species
            function listTheSpecies(fid,list)
                m = ceil(length(list)/4);
                r = length(list);
                n = 1;
                fprintf(fid,'!\n');
                for k = 1:m
                    fprintf(fid,'!');
                    for l = 1:4,
                        if n <= r
                            fprintf(fid,['\t',list{n}]);
                            n = n+1;
                        end
                    end
                    fprintf(fid,'\n');
                end
                fprintf(fid,'! \n');
            end

        % looking for efficencies  lines
            function b = isEfficienciesLine(cline)
                % an efficiencies line uses a slash as delimiter 
                b = ~isempty(strfind(cline,'/'));        
            end

        % remove a Species from the efficiencies line
            function cline = removeSpeciesFromEfficiencies(cline,deadSpeciesNamesList)
                for k = 1:length(deadSpeciesNamesList)
                    startSpecies = min(strfind(cline,[deadSpeciesNamesList{k},'/']));
                    if ~isempty(startSpecies)              
                        endspecies = min(strfind(cline(startSpecies:end),'/'));
                        endspecies = startSpecies+endspecies;
                        endspecies = endspecies+min(strfind(cline(endspecies+1:end),'/'));
                        cline(startSpecies:endspecies) = blanks(endspecies-startSpecies+1);
                    end       
                end
            end
        % fill with blanks
            function str = blank(n)
                    str = [];
                    for kk = 1:n
                        str = [str,' '];
                    end
            end
            function b =  efficienciesContainsDeadSpecies(cline,deadSpeciesNamesList)
                b = false;
                for k = 1:length(deadSpeciesNamesList)
                    lngthSpecName = length(deadSpeciesNamesList{k});
                    for l = 1:length(cline)-lngthSpecName-1;
                        foundSpec = strcmp(cline(l:l+lngthSpecName-1),deadSpeciesNamesList{k});
                        b = isLeft(cline,l) && isRight(cline,l) && foundSpec ; 
                        if b 
                            return
                        end               
                    end
                end 
                % local functions for identifying species by name
                function b = isLeft(cline,l)
                    % checks if a species string is a whole species name or only a part
                    % of a species name: left of the name...
                    b = (l==1) || (strcmp(cline(l-1),' ')) ;
                end
                % of a species name: right of the name can be a blank, "<", "+" "=" "(+", or "/"...
                function b = isRight(cline,l)
                    b = ( strcmp(cline(l+lngthSpecName),' ') ||  strcmp(cline(l+lngthSpecName),'<') ||...
                        strcmp(cline(l+lngthSpecName),'+')||  strcmp(cline(l+lngthSpecName),'=') ||...
                        strcmp(cline(l+lngthSpecName:l+lngthSpecName+1),'(+') ||...
                        strcmp(cline(l+lngthSpecName),'/')) ; 
                end
            end
        end

        function removeSpeciesCTI(mechFileIn,mechFileOut,deadSpeciesIndexList,deadSpeciesNamesList)
            % removeDeadSpecies(mechFileIn,mechFileOut,deadSpeciesIndexList)
            %   Removes the species (and Reactions) defined in deadSpeciesIndexList and
            %   deadSpeciesNamesList from mechFileIn, the resulting mechanism will be
            %   written in mechFileOut. 
            %   The input arguments deadSpeciesIndexList,deadSpeciesNamesList must be
            %   consistent, i.e. it holds  
            %          species(deadSpeciesIndexList) = deadSpeciesNamesList.
            
            
            
            % check for stupid inputs/output file declaration
            if strcmp(mechFileIn,mechFileOut)
                MException(...
                    'removeSpeciesCTI:SameInputOutputFile',...
                    'The input is the same as the output file. This is not allowed.').throw;
            end
            % try to open the files
            
            fprintf('removeSpeciesCTI: Indexing input file\n')
            try 
                indexCti = chemistry.indexCtiFile(mechFileIn);
            catch ME
                throw(ME)     
            end
            fprintf(['removeSpeciesCTI: Mechanism File contains ',...
                num2str(length(indexCti.species)),' species ',...
                num2str(length(indexCti.reactions)+...
                length(indexCti.three_bodies)+length(indexCti.falloffs)),...
                ' reactions.\n'])

            
            % check the input list
            if ~isnumeric(deadSpeciesIndexList)
                error('removeSpeciesCTI:ListMustBeNumeric',...
                    'Dead Species Index List must contain the INDICES of dead species')
            else
                % list is ok
                fprintf('Check input list: O.K.\n')
            end
            fprintf(['removeSpeciesCTI: There are ',num2str(length(deadSpeciesIndexList)),' species in dead species list.\n'])
            fprintf('removeSpeciesCTI: I found \n\t')
            listTheSpecies(deadSpeciesNamesList)
            fprintf(['removeSpeciesCTI: I will now try to delete all',...
                ' these species and related reacion.\n\t',...
                ' It will take only a few of your earth seconds.\n']);
            %----------------------------------------------------
            % some logical field to mark the "areas of action"
            
            % specify the lines we have to work or not to work
            % initialize all true
            nothingToDo = logical(1:indexCti.nlines);
            
            handleReactions = logical(1:indexCti.nlines); 
            handleReactions = ~handleReactions;
            handleIdealGas = logical(1:indexCti.nlines);
            handleIdealGas = ~handleIdealGas;
            notToCopy = logical(1:indexCti.nlines);
            notToCopy = ~notToCopy;




            % but we have to modify all ideal_gas blocks, mark the "to removed species"
            nothingToDo(indexCti.ideal_gas(1):indexCti.species(1)-2) = false;
            % we have an action for the to remove species..
            for k = 1:length(deadSpeciesIndexList)
                if deadSpeciesIndexList(k)==length(indexCti.species)
                    % the last species is choosen, then we mark all lines until the
                    % reactions block starts 
                    nothingToDo(indexCti.species(deadSpeciesIndexList(k)):...
                        min([indexCti.falloffs(1),indexCti.three_bodies(1),indexCti.reactions(1)])-8) = false;
                else    % the regular case
                    nothingToDo(indexCti.species(deadSpeciesIndexList(k)):indexCti.species(deadSpeciesIndexList(k)+1)-1) = false;
                end
            end
            
            % and for all reactions
            nothingToDo(indexCti.reactions) = false; % mark reactions be a canditate   
            nothingToDo(indexCti.reactions+1) = false; % mark the options line, it can be empty                                           
            nothingToDo(indexCti.three_bodies) = false; % mark reaction to  be a canditate 
            nothingToDo(indexCti.three_bodies+1) = false; % mark the efficiencies line                                                                                              
            nothingToDo(indexCti.falloffs) = false;       % mark falloffs to be a canditate
            nothingToDo(indexCti.falloffs+1) = false;     % there are four line after the falloff_reaction keyword
            nothingToDo(indexCti.falloffs+2) = false;     % we have to test
            nothingToDo(indexCti.falloffs+3) = false;
            nothingToDo(indexCti.falloffs+4) = false; 
            handleIdealGas(indexCti.ideal_gas(1):indexCti.species(1)-2) = true;
            % the special actions:
            
            % now we remove all species defined in  the parameter "deadSpeciesIndexList"
            % thew simplest action: to handle the "to remove species" 
            for k = 1:length(deadSpeciesIndexList)
                if deadSpeciesIndexList(k)==length(indexCti.species)
                    % the lst species is choosen, then we mark all lines until the
                    % reactions block starts
                    notToCopy(indexCti.species(deadSpeciesIndexList(k)):...
                        min([indexCti.falloffs(1),indexCti.three_bodies(1),indexCti.reactions(1)])-8) = true;
                else    % the regular case
                    notToCopy(indexCti.species(deadSpeciesIndexList(k)):indexCti.species(deadSpeciesIndexList(k)+1)-1) = true;
                end
            end
            handleReactions(indexCti.reactions) = true;
            handleReactions(indexCti.three_bodies) = true;
            handleReactions(indexCti.falloffs) = true;
            inSpecies = false;
            deadReactionRemove = false;
            threeBodyRemove = false;
            falloffRemove  = false;
            
            % open the files
            [fidMech,mesg] = fopen(mechFileIn,'r');
            if fidMech < 0 
                fprintf(mesg)
            end
            [fidNewMech,mesg] = fopen(mechFileOut,'w');
            if fidNewMech < 0 
                fprintf(mesg)
            end
            
            % for counting the actions
            nToRemoveReactions = 0;
            % a counter to make the ideal_gas mames unique
            nGasBlocks = 0;
            
            for k = 1:indexCti.nlines
                if nothingToDo(k)   
                    %-----------------------------------
                    % comments living species etc
                    % just copy...
                    %-----------------------------------
                    cline = fgetl(fidMech); 
                    fprintf(fidNewMech,[cline,'\n']);
                    %
                elseif handleIdealGas(k)
                    %-----------------------------------
                    % the ideal gas block contains also a species block that should be
                    % modified. We copy all lines exept this species block. However, its
                    % length is not fixed so we must scan it
                    % The species block starts with the keyword "species" and ends with the
                    % keyword "reactions". The keywords start at pos 7 in every line
                    cline = fgetl(fidMech); 
                    if length(cline)>=13 && strcmp(cline(7:13),'species')
                        inSpecies = true; % switch species block "on"
                    elseif length(cline)>=15&& strcmp(cline(7:15),'reactions')
                        inSpecies = false; % switch species block "off"
                    end
                    if inSpecies
                        % the "species action"
                        % analyse the line and remove the dead species
                        nline = chemistry.removeDeadSpeciesFromIdealGasBlock(cline,deadSpeciesNamesList);
                        fprintf(fidNewMech,[nline,'\n']);            
                    elseif length(cline)>=13 && strcmp(cline(7:14),'kinetics')
                        % can be omitted                %
                    elseif length(cline)>=13 && strcmp(cline(7:15),'transport')           % 
                        fprintf(fidNewMech,[cline,'\n']);   
                    elseif length(cline)>=13 && strcmp(cline(1:5),'ideal')
                        nGasBlocks = nGasBlocks + 1;
                        fprintf(fidNewMech,['ideal_gas(name = "gas_',num2str(nGasBlocks),'",\n']);              
                    else
                        fprintf(fidNewMech,[cline,'\n']);
                    end
                elseif notToCopy(k)
                    %-----------------------------------
                    % The  definition lines of dead species:
                    % removing by "comment out"
                    %-----------------------------------
                    %         cline = fgetl(fidMech);
                    fgetl(fidMech); 
                    %         fprintf(fidNewMech,['#  ',cline,'\n']);  
                elseif handleReactions(k)
                    %-----------------------------------
                    % The definition lines  of dead recations: 
                    % removing  by "comment out" 
                    %-----------------------------------
                    cline = fgetl(fidMech);
                    
                    % ckecking if the reaction contains a dead species
                    if chemistry.isDeadReaction(cline,deadSpeciesNamesList)
                        %             fprintf(fidNewMech,['#  ',cline,'\n']);
                        deadReactionRemove = true;
                        nToRemoveReactions = nToRemoveReactions +1;
                        %            
                        if strcmp(cline(1:5),'three') 
                            % mark the line after a dead three body as "to be removed"
                            threeBodyRemove = true;
                        elseif strcmp(cline(1:4),'fall') 
                            falloffRemove = true;
                        end
                    else
                        % if the reaction contains only living species, we copy it into
                        % the reduced mech file
                        deadReactionRemove = false;
                        falloffRemove = false;
                        threeBodyRemove = false;
                        fprintf(fidNewMech,[cline,'\n']);
                    end
                else
                    cline = fgetl(fidMech);
                    if length(cline)>16&&strcmp(cline(10:16),'options')&&deadReactionRemove
                        %             fprintf(fidNewMech,['#',cline,'\n']);  
                        deadReactionRemove = false;
                    elseif threeBodyRemove
                        % remove the efficnecies line of a three body reaction
                        threeBodyRemove = false;
                        %             fprintf(fidNewMech,['#',cline,'\n']);           
                    elseif falloffRemove
                        % remove the falloff lines 2-5
                        % analyse the first two characters in the lines
                        if length(cline)>=11
                            key = cline(10:11);
                        else
                            key = ' ';
                        end
                        switch key
                            case {'kf' 'fa'}
                                %                     fprintf(fidNewMech,['#',cline,'\n']);                    
                            case 'ef'
                                %                     fprintf(fidNewMech,['#',cline,'\n']);   
                                falloffRemove = false;
                            otherwise
                        end
                    elseif length(cline)>=21,
                        if strcmp(cline(10:21),'efficiencies') 
                            % remove the dead species from the efficiencies list
                            newEfficiencies = chemistry.removeDeadFromEfficiencies(cline,deadSpeciesNamesList);
                            fprintf(fidNewMech,[newEfficiencies,'\n']);
                        else
                            fprintf(fidNewMech,[cline,'\n']);
                        end
                    else
                        fprintf(fidNewMech,[cline,'\n']);
                    end
                    % up to now we write a empty line instead to make the "real action"
                    % here we define the action    
                end
            end
            fprintf('removeDeadSpecies: Try to write the files and close them.\n')
            % close all files
            
            st = fclose(fidNewMech);
            if st<0
                MException('removespeciesCTI:fcloseFails',...
                    'fclose fails to close the mechfile').throw;
            end
            st = fclose(fidMech);
            if st<0
                MException('removespeciesCTI:fcloseFails',...
                    'fclose fails to close the mechfile').throw;
            end
            st = fclose('all');
            if st<0
                MException('removespeciesCTI:fcloseFails',...
                    'fclose fals to close the mechfile').throw;
            end
            
            fprintf(['removeDeadSpecies: Afterall, I remove ',...
                num2str(nToRemoveReactions),...
                ' of ',num2str(length(indexCti.reactions)+...
                length(indexCti.three_bodies)+length(indexCti.falloffs)),...
                ' reactions'])
            fprintf(['...or round ',...
                num2str(round(100*nToRemoveReactions/...
                (length(indexCti.reactions)+length(indexCti.three_bodies)...
                +length(indexCti.falloffs)))),' percent. \n'])
            fprintf('removeDeadSpecies: OK. Thats all, see you!\n\n')  
            
            

            % a tiny helper function: it lists a cell array full of strings
            function listTheSpecies(list)
                m = ceil(length(list)/4);
                r = length(list);
                n = 1;
                fprintf('\n')
                for kk = 1:m,
                    for l = 1:4,
                        if n <= r
                            fprintf(['\t',list{n}]);
                            n = n+1;
                        end
                    end
                    fprintf('\n')
                end
                fprintf('\n')
            end
        end
                 
        function species = scanSpeciesList(elementsLine)
            % utility function for removeSpecies
            % it reads the species from the elementsLine and gives back a cell array
            % containing the species names
            
            species = [];
            l = 0;
            candidate = [];
            spec = [];
            for k = 1:length(elementsLine)
                if  ~strcmp(elementsLine(k),' ')
                    spec = [spec,elementsLine(k)];
                end
                if strcmp(elementsLine(k),' ')
                    l = l+1;
                    if ~isempty(spec)
                        candidate{l}= spec;
                        spec = [];
                    end
                end
            end
            l = 1;
            for k = 1:length(candidate)
                if ~(isempty(candidate{k})||strcmp(candidate{k},'""",'))
                    species{l}=candidate{k};
                    l = l+1;
                end
            end            
        end
         
        function species = selectSpeciesFromList(elementsLine)
            % utility function for removeSpecies
            % it reads the species from the elementsLine and gives back a cell array
            % containing the species names
            l = 0;
            inList = false;
            spec = [];
            for k = 1:length(elementsLine)
                if inList && ~(strcmp(elementsLine(k),' ')||strcmp(elementsLine(k),'"'))
                    spec = [spec,elementsLine(k)];
                end
                if inList && (strcmp(elementsLine(k),' ')||strcmp(elementsLine(k),'"'))
                    l = l+1;
                    if ~isempty(spec)
                        candidate{l}= spec;
                        spec = [];
                    end
                end
                if strcmp(elementsLine(k),'"')
                    inList = swtch(inList);
                end
            end
            l = 1;
            for k = 1:length(candidate)
                if ~(isempty(candidate{k})||strcmp(candidate{k},'=>')||strcmp(candidate{k},'<=>')||strcmp(candidate{k},'+'))
                    species{l}=candidate{k};
                    l = l+1;
                end
            end
            function b = swtch(b)
                b = ~b;                
            end
        end
        
        function fractionsString = makeFractions(species,fractions)
        % fractionsString = makeMoleFractions(species,fractions)
        % Creates a mole or mass fractions string for the solution object
        % e.g. for use with the chemistry object
        % species is a cell array of species names and fractions is a vector of
        % length n of species. 
        % The fractions were normalized by fractions = fractions/sum(fractions).
        % (c) 2011 by U. Pruefert for Virtuhcon

            nspe = length(species);
            fractions = fractions/sum(fractions);
            fractionsString = [];
            for k = 1:nspe-1
                fractionsString = [fractionsString,species{k},':',num2str(fractions(k)),','];
            end
            fractionsString = [fractionsString,species{nspe},':',num2str(fractions(nspe))];
        end
        
        function moles = getMoleFractionsFromInitialFractions(str,varargin)
            % extracts molefractions from the (user input) inital mole
            % fractions and normalze them (if necessary)
            speciesMole =[];

            k = 1; 
            in = false;
            for i=1:length(str)   
                if strcmp(str(i),':')
                    in = true;               
                end
                if isnumber(str(i)) && in
                    speciesMole = [speciesMole,str(i)] ; 
                end
                if (strcmp(str(i),',')||i==length(str)||strcmp(str(i),' ')) && in
                    moles(k) = str2double(speciesMole);
                    k=k+1;
                    in = false;
                    speciesMole=[];
                end
            end

            switch nargin
                case 3
                    if strcmp(varargin{1},'normalize')
                        if strcmp(varargin{2},'off')  
                            %                
                        elseif strcmp(varargin{2},'on')
                            moles = moles./sum(moles);
                        end
                    end
                otherwise
                    moles = moles./sum(moles);
            end






            function b = isnumber(char)
                % checks if char is a number, i.e. it contains '0'...'9' and '.'
                b = ((double(char)<58 && double(char)>47))||double(char)==46||double(char)==45||double(char)==101;
            end
        end
        
        
    end
        
    methods(Abstract)
        obj2 = copy(obj1);
    end    
    
    methods(Static,Hidden)
        addlistener
        eq 
        findprop 
        gt   
        le   
        ne       
        % delete  
        findobj
        ge  
        % isvalid 
        lt     
        notify       
    end
end