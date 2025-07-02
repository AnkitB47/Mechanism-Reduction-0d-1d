classdef nasa
    % nasa polynomial class 
    % (C)  U. Pruefert for VIRTUHCON 2011
    % Uses XML4MAT by Jonas S. Almeida, Shuyuan Wu, and Eberhard O. Voit for
    % evaluating the cantera xml format.
    
    properties(SetAccess=private)
        speciesName;
        coeff;
        T0 = 298.15; % standard-state (SATP: Standard Ambient Temperature
                     % and Pressure) T0 = 298.15 K
    end
    
    methods
        function obj = nasa(varargin)
            % contructor for nasa class
            % Usuage:
            % obj = nasa(xmltree)
            % obj = nasa(path_to_cxml_file) 
            % obj = nasa(mat-file)
            % The input MUST be a valid cxml file or a xml tree that comes
            % from such a file. xml tree can be a variable in the workspace
            % or stored in a mat file.
            % (C)  U. Pruefert for VIRTUCHON 2011
            switch length(varargin)
                case 0
                    obj.speciesName = [];
                    obj.coeff = [];             
                case 1
                    %read the mech xml-file and write it into a structure
                    if  isstruct(varargin{1})
                        xmlTree = varargin{1};
                        if isempty(xmlTree(1).speciesData)
                            error('speciesData:NoXMLTree','The xml Tree input is incorrect.')
                        end
                    else
                        try 
                            load(varargin{1});
                        catch ME
                            fprintf([ME.message,'\n']);
                            try
                                if ~exist( 'xmlTree','var')
                                    xmlTree = xml2struct(varargin{1});
                                end
                            catch ME
                                error(ME.identifier,ME.message)
                            end
                        end
                    end
                    % the number of species
                    nSpecies = length(xmlTree(1).speciesData);
                    for k = 1:nSpecies
                        %  the name
                        obj.speciesName{k} = xmlTree(1).speciesData(k).species.name;
                        % the number of NASA Polynomes, 1..2
                        nNasa = length(xmlTree(1).speciesData(k).species.thermo);
                        switch nNasa
                            case 1
                                obj.coeff(k).splitt = inf;
                                % the coefficent are strings so we convert
                                % them into numerics. str2 num take care for the fact that it is a vector... 
                                obj.coeff(k).low = str2num(xmlTree(1).speciesData(k).species.thermo(l).NASA.floatArray.floatArray);
                            case 2 
                                obj.coeff(k).low = str2num(xmlTree(1).speciesData(k).species.thermo(1).NASA.floatArray.floatArray); %#ok<*ST2NM>
                                obj.coeff(k).high = str2num(xmlTree(1).speciesData(k).species.thermo(2).NASA.floatArray.floatArray);
                                obj.coeff(k).splitt = xmlTree(1).speciesData(k).species.thermo(1).NASA.Tmax;                            
                            otherwise
                                % fuck!!!
                        end
                    end
                otherwise
            end
        end
        % -------  nasaValue method -------
        function [varargout] = nasaValue(obj,varargin)
            % Evaluates the nasa polynomial at  temperature T
            % calls can be 
            % cp = nasa.nasaValue(T)  all heat capacities cp (row vector) 
            % cp = nasa.nasaValue(speciesIndex, T) heat capacities cp for
            % the species speciesIndex
            % cp = nasa.nasaValue(speciesName,T) heat capacities cp for the
            % species speciesName (String) If speciesName is unkown, or
            % speciesIndex is larger than the number of species,
            % nasaValue aborts with an error.
            % Outputs can be the heat capacity only, heat capacity and
            % Enthalpie and heat capacity, Enthalpy and Entropy:
            % cp = nasa.nasaValue(T)   
            % [cp,h] = nasa.nasaValue(T)
            % [cp,h,s] = nasa.nasaValue(T)
            %
            % Well, if your using Matlab >2009 you can call
            % [~,h,~] = nasa.nasaValue(T) to obtain enthalpy only ...
            % or  [~,~,s] = nasa.nasaValue(T)  for entropy only
            % 
            % (C)  U. Pruefert for VIRTUCHON 2011
            
            switch nargin
                case 2
                   species = 1:length(obj.speciesName);
                   temp = varargin{1};
                case 3
                    species = varargin{1}; 
                    if isnumeric(species)
                        species = double(species);
                    end
                    temp = varargin{2};
                otherwise
                    error('nasaValue:WrongNumberArguments','Wrong number of input arguments, try help nasaValue.')
            end
            %
            switch class(species)
                case 'double'
                    cp = zeros(1,length(species));
                    h0 = zeros(1,length(species));
                    s0 = zeros(1,length(species));
                    for k = 1:length(species)
                        temp2 = str2double(obj.coeff(species(k)).splitt);
                        if temp < temp2
                            % low Temp
                            cp(k) = obj.coeff(species(k)).low(1:5)*temp.^(0:4)'; % heat capacity
                            h0(k) = 1./(1:5).*obj.coeff(species(k)).low(1:5)*temp.^(0:4)'+obj.coeff(species(k)).low(6)/temp; % enthalpie
                            s0(k) = obj.coeff(species(k)).low(1)*log(temp)+obj.coeff(species(k)).low(7)+... % entropie
                               1./(1:4).*obj.coeff(species(k)).low(2:5)*temp.^(1:4)';
                        else
                            % high Temp                    
                            cp(k) = obj.coeff(species(k)).high(1:5)*temp.^(0:4)';
                            h0(k) = 1./(1:5).*obj.coeff(species(k)).high(1:5)*temp.^(0:4)'+obj.coeff(species(k)).high(6)/temp;
                            s0(k) = obj.coeff(species(k)).high(1)*log(temp)+obj.coeff(species(k)).high(7)+...
                               1./(1:4).*obj.coeff(species(k)).high(2:5)*temp.^(1:4)';
                        end
                    end
                case 'char'
                    r = 0;
                    for k = 1:length(obj.speciesName)                                          
                        if strcmp(obj.speciesName{k},species)                                                  
                            break; 
                        else
                            r = k+1;
                        end
                    end
                    %
                    if r>length(obj.speciesName)
                            error('nasaValue:SpeciesNotFound',['Fatal: Species ', species,' not found'])
                    end
                    temp2 = str2double(obj.coeff(k).splitt);
                    if temp < temp2
                        % low Temp
                        arr = obj.coeff(1,k).low(1:5); % TODO
                        cp = arr*temp.^(0:4)';
                        h0 = 1./(1:5).*obj.coeff(k).low(1:5)*temp.^(0:4)'+obj.coeff(k).low(6)/temp;
                        s0 = obj.coeff(k).low(1)*log(temp)+obj.coeff(k).low(7)+...
                           1./(1:4).*obj.coeff(k).low(2:5)*temp.^(1:4)';
                    else
                        % high
                        cp = obj.coeff(k).high(1:5)*temp.^(0:4)';
                        h0 = 1./(1:5).*obj.coeff(k).high(1:5)*temp.^(0:4)'+obj.coeff(k).high(6)/temp;
                        s0 = obj.coeff(k).high(1)*log(temp)+obj.coeff(k).high(7)+...
                           1./(1:4).*obj.coeff(k).high(2:5)*temp.^(1:4)';
                    end
                otherwise
                % fuck again
                error('nasaValue:WrongInputFormat','The species must be given as string or double (vector).')
           end
           switch nargout
               case 0
                   %
               case 1
                   varargout{1} = cp;
               case 2
                   varargout{1} = cp;
                   varargout{2} = h0;
               case 3
                   varargout{1} = cp;
                   varargout{2} = h0;
                   varargout{3} = s0;
               otherwise
                   % fuck the third
                   error('nasaValue:WrongNumberArguments','Wrong number of output arguments, try help nasaValue.')
           end
        end

% ------ nasaValueIntegral method  -------
        function [varargout] = nasaValueIntegral(obj,varargin)
            % Evaluates the integral of the nasa polynomial for the molar
            % heat capacity temperature T, i.e. :
            % 
            % int_{T_0}^T c_p(X) dX
            % calls can be 
            % cpInt = nasa.nasaValueIntegral(T) all integrals for the heat capacities cp (row vector) 
            % cpInt = nasa.nasaValueIntegral(speciesIndex, T) integral of the heat capacity cp for the species speciesIndex
            % cpInt = nasa.nasaValueIntegral(speciesName,T) integral of the heat capacity cp for the
            % species speciesName (String) If speciesName is unkown, or
            % speciesIndex is larger than the number of species,
            % nasaValueIntegral aborts with an error
            % (C)  Philipp Maeser for VIRTUCHON 2011
            switch nargin
                case 2
                   species = 1:length(obj.speciesName);
                   temp = varargin{1};
                case 3
                    species = varargin{1}; 
                    if isnumeric(species)
                        species = double(species);
                    end
                    temp = varargin{2};
                otherwise
                    error('nasaValueIntegral:WrongNumberArguments','Wrong number of input arguments, try help nasaValueIntegral.')
            end
            %
            switch class(species)
                case 'double'
                    cpInt = zeros(1,length(species));
                    for k = 1:length(species)
                        temp2 = str2double(obj.coeff(1,species(k)).splitt);
                        if temp < temp2
                            cpInt(k) = 0;
                            % low temperature
                            for i=1:5
                                cpInt(k) = cpInt(k) + (obj.coeff(1,species(k)).low(i)*(temp^i-obj.T0^i))/i;
                            end
                        else
                            % high temperature              
                            for i=1:5
                                cpInt(k) = cpInt(k) + (obj.coeff(1,species(k)).high(i)*(temp^i-obj.T0^i))/i;
                            end
                        end
                    end
                case 'char'
                    r = 0;
                    for k = 1:length(obj.speciesName)                                          
                        if strcmp(obj.speciesName{k},species)                                                  
                            break; 
                        else
                            r = k+1;
                        end
                    end
                    %
                    if r>length(obj.speciesName)
                            error('nasaValueIntegral:SpeciesNotFound',['Fatal: Species ', species,' not found'])
                    end
                    cpInt = 0;
                    temp2 = str2double(obj.coeff(k).splitt);
                    if temp < temp2
                        % low temperature
                        for i=1:5
                            cpInt = cpInt + (obj.coeff(1,k).low(i)*(temp^i-obj.T0^i))/i;
                        end
                    else
                        % high temperature
                        for i=1:5
                            cpInt = cpInt + (obj.coeff(1,k).high(i)*(temp^i-obj.T0^i))/i;
                        end
                    end
                otherwise
                % fuck again
                error('nasaValueIntegral:WrongInputFormat','The species must be given as string or double (vector).')
           end
           switch nargout
               case 0
                   %
               case 1
                   varargout{1} = cpInt;
               otherwise
                   % fuck the third
                   error('nasaValueIntegral:WrongNumberArguments','Wrong number of output arguments, try help nasaValueIntegral.')
           end
        end
    end
end