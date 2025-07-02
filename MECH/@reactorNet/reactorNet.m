classdef reactorNet
% reactorNet class for simulation a reactor net.           
% Two reserviors and two reactors are coupled by five mass flow controllers.    
% U.Pruefert 2011
    
    
    
    properties (SetAccess = private)
       % 
    end
    
    properties (SetAccess = public)        
        % only for reactor networks needed: Mass flow controlers
        reactor1;
        reactor2;
        reservoir1;
        reservoir2;
        massFlowControllers;        
    end
    
      
    
    methods
        % The CONSTRUCTOR
        function obj = reactorNet(varargin)
            switch nargin 
                case 0
                    mechfile = [];
                case 1
                    mechfile = varargin{1};
                otherwise
            end
            % two batchReactor objects
            obj.reactor1 = batchReactor(mechfile);
            obj.reactor2 = batchReactor(mechfile);
            
            % two reservoirs. reservoir is here only a container for the
            % initial conditions. We do not need all of the batchReactor
            % stuff (in contrast to the Cantera implementation)
            obj.reservoir1 = batchReactor(mechfile);
            obj.reservoir2 = batchReactor(mechfile);          
                       
            % only for reactor networks needed:
            % it seems to be enough to have the values for the mass flow
            % stored (?)
            obj.massFlowControllers = zeros(1,5);            
        end
        
        % The DESTRUCTOR
        function obj = clear(obj)
            try                 
                cleanup                
            catch ME
                fprintf([ME.message,'\n'])
            end
        end
        
        % a plot function
        function plotMoleFractions(obj,varargin)
            subplot(1,2,1)
                switch length(varargin)
                    case 0
                        obj.reactor1.plotMoleFractions;
                    case 1
                        obj.reactor1.plotMoleFractions(varargin{1});
                    case 2
                        obj.reactor1.plotMoleFractions(varargin{1},varargin{2});
                    otherwise
                end    
            subplot(1,2,2)
                switch length(varargin)
                    case 0
                        obj.reactor2.plotMoleFractions;
                    case 1
                        obj.reactor2.plotMoleFractions(varargin{1});
                    case 2
                        obj.reactor2.plotMoleFractions(varargin{1},varargin{2});
                    otherwise
                end   
        end
        
        %%%%%%%%%% methods block ends here %%%%%%%%%%
    end
end