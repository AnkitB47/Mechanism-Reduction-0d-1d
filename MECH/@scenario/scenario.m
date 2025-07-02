classdef scenario
    % 
    
    properties(SetAccess=private)
        temperatures;
        pressures;
        moleFractions;
    end
    
    methods
        function obj = scenario
            obj.temperatures = [];
            obj.pressures = [];
            obj.moleFractions = [];
        end
        
        function obj = setTemperatures(obj,t0,t1,n)
            switch nargin
                case 2
                    if length(t0)<2
                        fprintf('Set a single temperature.')
                    end
                    obj.temperatures = t0;             
                case 3                    
                    obj.temperatures = linspace(t0,t1,3);
                case 4                    
                    obj.temperatures = linspace(t0,t1,n);
                otherwise
                    %
                    error('setTemperatures:WrongNumberInputArguments','Wrong number of input arguments. It must be 1,2 or 3.')
            end
             
           
        end
        
        function obj = setPressures(obj,p0,p1,n)
            switch nargin
                case 2
                    if length(p0)<2
                       fprintf('Set a single pressure.\n') 
                    end
                    obj.pressures = p0;                     
                case 3                    
                    obj.pressures = linspace(p0,p1,3);
                case 4                    
                    obj.pressures = linspace(p0,p1,n);
                otherwise
                    %
                    error('setsetPressures:WrongNumberInputArguments','Wrong number of input arguments. It must be 1,2 or 3.')
            end
        end
        
        function obj = addMoleFractions(obj,moleFraction)            
            obj.moleFractions{length(obj.moleFractions)+1} = moleFraction;  
        end
        
        function obj = delMoleFractions(obj,toDelete)
            l = 0;
            for k = 1:length(obj.moleFractions)
                if ~(k==toDelete)
                    l = l+1;
                    intermediate{l} = obj.moleFractions{k};
                end
            end
            obj.moleFractions = intermediate;
        end
    end
    
end

