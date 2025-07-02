classdef batchReactor < chemistry
    % BATCHREACTOR class.
    % contructor call: batchReactor(mechfile) 
    % U.Pruefert 2011
    
    
    % private = results of the simulation
    properties (SetAccess = private)
        % only 0D-reactor specific properties are defined here
        % Note that mass and mole fractions are come from the
        % chemistry super-class 
        
        times;  % timesteps
        lambda; % times scales
        %-----------------------------
        ignt;
        %-----------------------------
        
    end
    
    % public = user data 
    properties (Access = public)          
        % type of reactor: constant volume, constant pressure,
        % common reactor, the tibox-reactor, etc.    
        reactor;
        %-----------------------------
        initialMoleFractions;
        initialMassFractions; 
        %-----------------------------
        initialTemperature;
        initialPressure;
        %-----------------------------
        simulationTime;
        equilibrium = 'false';
        %-----------------------------
        % for the common reactor...
        % if not set, standard values are used
        heatflux;
        vdot;
        area;
        %-----------------------------        
    end
    
  
    
    methods(Access=public)
        % The CONSTRUCTOR
        function obj = batchReactor(varargin)
        if nargin == 0 
            % Create empty cell array if no input argsuments
            super_class_args = {}; 
        else
            % Use specified arguments
            super_class_args{1} = varargin{1};           
        end
        % batchReactor(mechanism)         
        obj = obj@chemistry(super_class_args{:});  
        
        end 
        
        % ...some signatures of functions  
        solve(obj,varargin); 
        setInitialMassFractions(obj,names,values);  
        speciesList = getInitialSpecies(obj);
        varargout = switchReactions(obj,reactions); 
        varargout = timeScaleJacobi(obj,flag1,val1);
        varargout = timeScale(obj,varargin)
        varargout = speciesProductionRates(obj,varargin);
        obj = copy(obj); % needs argument!
        M = sensitivityMatrix(obj,y);
        varargout = equilibrate(obj,arg);
    end
    
    methods(Access=public,Hidden)
        % time scales for timeScale method...
        varargout = ijts(obj,flag1,val1);
        varargout = evts(obj,flag1,val1,flag2,val2);
        varargout = timeScaleCaudalMod(obj,flag1,val1);
        varargout = pvts(obj,species);
        varargout = gpvts(obj,species);        
    end
    
    methods (Access=private,Hidden)
        % source terms for batch reactors...
        w = conhp(obj,t,y);
        w = conuv(obj,t,y);
        w = conht(obj,t,y);
        w = creactor(obj,t,y);  
        
        % reactor types - called by solve...
        commonReactor(obj);
        constTempReactor(obj);
        constPressReactor(obj);
        constPressNetReactor(obj);
        constVolumeReactor(obj);
        tibox(obj)  
        % helper for time scale computing
        [df,f] = sensitivitY(obj,y);        
        setTimes(obj,times);      
        setTemp(obj,times);      
        setDens(obj,dens);      
        setPres(obj,pres);
    end
    
    
end
