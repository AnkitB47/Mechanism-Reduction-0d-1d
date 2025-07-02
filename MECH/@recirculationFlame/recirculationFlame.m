classdef recirculationFlame
    properties (SetAccess = public)
        reactor;
        flame;
        ratio;
        recirculatingMixture;        
    end
    
    methods
        function obj = recirculationFlame(mech,mixModel)
            %constructor method
            % obj = recirculationFlame(mech) 
            % obj = recirculationFlame(mech,mixModel)
            switch nargin
                case 0
                    error('recirculationFlame:recirculationFlame:numberInputs',...
                        'Creating recirculationFlame Object failed: No mechanism given.')
                case 1
                    obj.reactor = batchRector(mech);
                    obj.flame = flame1d(mech);
                case 2
                     obj.reactor = batchReactor(mech);
                     obj. flame = flame1d(mech,mixModel);
                otherwise
                    error('recirculationFlame:recirculationFlame:numberInputs',...
                        'Wrong number of inputs.')                    
            end
            % this flame contains a counterflow flame and a constant
            % pressure reactor
            obj.flame.flameType = 'counterFlowFlame';
            obj.reactor.reactor = 'constPressReactor' ; 
        end
        
        function obj = solve(obj)
            % solves the recirculation Model
            % (i)   solve the counterflow flame
            % (ii)  measure the mass fractions in the flame zone
            % (iii) mix with the recirculation streme
            % (iv) solve the reactor
            
            
            fprintf('########################################\n')
            tic
            obj.flame.solve; 
            neededtime = toc;
            fprintf(['#\t\t\t flame solved in ',num2str(neededtime),' sec\t\t\t\n#\n'])
            fprintf('########################################\n')
            Y = obj.flame.massFractions(:,obj.flame.locateFlameFront);
            obj.reactor.initialMoleFractions = ...
                obj.reactor.mixFractions(obj.reactor.vector2CanteraString(obj.flame.getLivingSpeciesNames,...
                                                                        Y(obj.flame.getLivingSpeciesIndex)),...
                                        obj.reactor.mole2mass(obj.recirculatingMixture),obj.ratio); % in MOLE FRACTIONS !!!
            if isempty(obj.reactor.initialTemperature)
%                 obj.reactor.initialTemperature = obj.flame.meanTemperature;      
                obj.reactor.initialTemperature =  obj.mixtureTemperature(obj.reactor.vector2CanteraString(obj.flame.getLivingSpeciesNames,...
                                                                        Y(obj.flame.getLivingSpeciesIndex)),...
                                        obj.reactor.mole2mass(obj.recirculatingMixture),max(obj.flame.temp),obj.flame.tInLeft,obj.ratio);
            end
            tic      
            obj.reactor = obj.reactor.solve;   
            neededtime = toc; 
            fprintf(['#\t\t\t reactor solved in ',num2str(neededtime),' sec\t\t\t\n#\n'])
            fprintf('########################################')
        end
        
        function plotMoleFractions(obj,varargin)
            switch length(varargin)
                case 1
                    subplot(2,1,1)
                    obj.flame.plotMoleFractions;
                    subplot(2,1,2)
                    obj.reactor.plotMoleFractions;
                case 2
                    subplot(2,1,1)
                    obj.flame.plotMoleFractions(varargin{1});
                    subplot(2,1,2)
                    obj.reactor.plotMoleFractions(varargin{1});
                case 3
                    subplot(2,1,1)
                    obj.flame.plotMoleFractions(varargin{1},varargin{2});
                    subplot(2,1,2)
                    obj.reactor.plotMoleFractions(varargin{1},varargin{2});
                otherwise
            end
        end
        
        function plotMassFractions(obj,varargin)
            switch length(varargin)
                case 0
                    subplot(2,1,1)
                    obj.flame.plotMassFractions;
                    subplot(2,1,2)
                    obj.reactor.plotMassFractions;
                case 1
                    subplot(2,1,1)
                    obj.flame.plotMassFractions(varargin{1});
                    subplot(2,1,2)
                    obj.reactor.plotMassFractions(varargin{1});
                case 2
                    subplot(2,1,1)
                    obj.flame.plotMassFractions(varargin{1},varargin{2});
                    subplot(2,1,2)
                    obj.reactor.plotMassFractions(varargin{1},varargin{2});
                otherwise
            end
        end
        
        
        function plotTemp(obj)
            subplot(2,1,1)
            obj.flame.plotTemp;
            subplot(2,1,2)
            obj.reactor.plotTemp;
        end
    end  
    methods(Static)
    end
end