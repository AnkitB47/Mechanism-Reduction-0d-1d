classdef speciesData
    %speciesData class
    % contains all species related data given by the mechanism 
    % speciesData(XML-FilE)
    % speciesData(XMLTREE)
    % Properties
    %     molecularWeights (g/mol)
    %     heatCapacity(T)  (dimension less)
    %     enthalpy(T)      (dimension less)    
    % Some tests:
    % molecularWeight gives the same result as molecularWeights(gas) with
    % an error less of 1.0e-13 (molecular weight is in kg/mol)
    % Error in the Enthalpies [~,h,~] = obj.ethalpy(300) less
    % than 1.0e-12. (comp. with the related CANTERA call)
    % 
    % (c) Uwe Pruefert for VIRTUCHON 2011.
    
    properties
        % public
    end
    
    properties(SetAccess=private)        
        molecularWeights;       
    end
    
    properties(Access=private)
        nasaData;        
    end
    
    methods
        function obj = speciesData(varargin)
            % constructor
            atomMasses{1}.element = 'H';
            atomMasses{1}.value = 1.00794;
            
            atomMasses{2}.element  = 'O';
            atomMasses{2}.value = 15.9994;
            
            atomMasses{3}.element  = 'N';
            atomMasses{3}.value = 14.00674;
            
            atomMasses{4}.element  = 'C';
            atomMasses{4}.value = 12.011;
           
            atomMasses{5}.element  = 'AR';
            atomMasses{5}.value =  39.948;
            
            atomMasses{6}.element  = 'HE';
            atomMasses{6}.value = 4.00260;
            
            atomMasses{7}.element  = 'NE';
            atomMasses{7}.value = 20.1797;
            
            
            switch nargin
                case 0
                case 1
                    if  isstruct(varargin{1})
                        xmlTree = varargin{1};
                        if isempty(xmlTree(1).speciesData)
                            error('speciesData:NoXMLTree','The xml tree input is incorrect.')
                        end
                    elseif ischar(varargin{1})
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
                    else
                        error('speciesData:WrongInputFormat','The input must be a xml file or an xml tree');
                    end
                    obj.molecularWeights = computeMolecularWeights(xmlTree);
                    obj.nasaData = nasa(xmlTree);
                otherwise
                    error('speciesData:WrongNumberArguments','The number of input arguments is wrong, try help speciesData');                    
            end
            
            function masses = computeMolecularWeights(xmlTree)
                %  
                % the molecular weight in kg/mol
                nSpecies = length(xmlTree(1).speciesData);
                masses = zeros(1,nSpecies);
                for k = 1:nSpecies
                    moleculeDefinitionString = xmlTree(1).speciesData(k).species.atomArray;
                    blankPos = strfind(moleculeDefinitionString,' ');
                    nElementsInMolecule = length(blankPos);
                    elementDefinitionString = [];
                    elementDefinitionString{1} = moleculeDefinitionString(1:blankPos(1)-1);  
                    for l = 2:nElementsInMolecule                                
                        elementDefinitionString{l} = moleculeDefinitionString(blankPos(l-1)+1:blankPos(l)-1);                            
                    end                   
                    for l = 1:length(elementDefinitionString)
                        splitter = strfind(elementDefinitionString{l},':');
                        element = elementDefinitionString{l}(1:splitter-1);
                        atoms = str2num(elementDefinitionString{l}(1+splitter:end)); %#ok<*ST2NM>
                        switch upper(element)
                            case 'H'
                                masses(k) = masses(k)+atomMasses{1}.value*atoms;
                            case 'O'
                                masses(k) = masses(k)+atomMasses{2}.value*atoms;
                            case 'N'
                                masses(k) = masses(k)+atomMasses{3}.value*atoms;
                            case 'C'
                                masses(k) = masses(k)+atomMasses{4}.value*atoms;
                            case 'AR'
                                masses(k) = masses(k)+atomMasses{5}.value*atoms;
                            case 'HE '
                                masses(k) = masses(k)+atomMasses{6}.value*atoms;  
                            case 'NE'
                                masses(k) = masses(k)+atomMasses{7}.value*atoms;  
                            otherwise
                                % fuck - fuck - fuck
                                error('speciesData:UnkownElement',['The element ',upper(element),' is inknown'])
                        end
                    end
                end
            end
        end
        
        function val = heatCapacity(obj,varargin)
            % val = heatCapacity(Temp)
            % Heat Capacity cp will be returend
            % species.nasaPolynom(Temperature)
            % species.nasaPolynom(SpeciesNumber,Temperature)
            % species.nasaPolynom(SpeciesName,Temperature)
            %
            % this function is a litle bit redundancy, the call cp = obj.heatCapacity(T)
            % is the same as cp = obj.nasaData.nasaValue(T)
            % Note: this is a dimesionless value, multiply it by R to
            % obtain it in J/kg/K etc.
            % (c) U.P. for VIRTUHCON 2011
            % Error lesser than 1.0.e-13 (comp. with the related CANTERA
            % call)            
            switch nargin
                case 1
                    % no temperature, set it to 293 K
                    val = obj.nasaData.nasaValue(293);
                case 2
                    val = obj.nasaData.nasaValue(varargin{1});
                case 3
                    val = obj.nasaData.nasaValue(varargin{1},varargin{2});
                otherwise
                    error('speciesData:WrongNumberOfArguments','Wrong number of input arguments');
            end           
        end
        %
        function val = enthalpy(obj,varargin)
            % val = enthalpy(Temp)
            % val = enthalpy(species,T)
            % The enthalpy h will be returend            
            % Note, that T is always th LAST argument.
            % This function is a litle bit redundant, the call cp = obj.enthalpy(T)
            % is the same as cp = obj.nasaData.nasaValue(T)
            % Note that this is "Pure species non-dimensional enthalpy". To
            % obtain a value in J/mol etc. multiply it by the Temperature
            % and the Gas-Constant R.
            % (c) U.P. for VIRTUHCON 2011
            % Error lesser than 1.0.e-13 (comp. with the related CANTERA
            % call)
            switch nargin
                case 1
                    % no temperature, set it to 293 K
                    val = obj.nasaData.nasaValue(293);
                case 2
                    [~,val] = obj.nasaData.nasaValue(varargin{1});
                case 3
                    [~,val] = obj.nasaData.nasaValue(varargin{1},varargin{2});
                otherwise
            end  
        end
        
        function val = entropy(obj,varargin)
            % val = entropy(Temp)
            % The etropy will be returend            
            % Note, that T is always th LAST argument.
           
            % Note that this is "Pure species non-dimensional entropy". To
            % obtain a value in J/mol etc. multiply it by the Temperature
            % and the Gas-Constant R.
            % (c) U.P. for VIRTUHCON 2011
            % Error lesser than 1.0.e-13 (comp. with the related CANTERA
            % call)
            switch nargin
                case 1
                    % no temperature, set it to 293 K
                    val = obj.nasaData.nasaValue(293);
                case 2
                    [~,~,val] = obj.nasaData.nasaValue(varargin{1});
                case 3
                    [~,~,val] = obj.nasaData.nasaValue(varargin{1},varargin{2});
                otherwise
            end  
        end     
        
       
        %
    end % methods   
    %
end

