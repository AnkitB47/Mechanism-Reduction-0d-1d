function obj = solveForIgnition(obj,varargin)
% obj = obj.solveForIgnition(marker) solves the reactor
% but w/o the postprocesing of commonReactor.
% Only times and temp were evaluated and the ignition time is computed.
% The marker should be a species which marks the ignition, like e.g. OH
obj.reactor = 'commonReactor';
obj = obj.clearFractions;
if  isa(obj.mechanism,'char') && isempty(obj.gas)
    obj.gas = IdealGasMix(obj.mechanism)  ;
    fprintf('chemistry.constVolumeReactor: Make obj.gas object\n')
elseif ~isIdealGas(obj.gas)   
    error('obj.gas is not an Ideal Gas Mix')
end

if isempty(obj.initialMoleFractions)
    set(obj.gas,'T',obj.initialTemperature,'P',obj.initialPressure,'Y',obj.initialMassFractions);
else
    set(obj.gas,'T',obj.initialTemperature,'P',obj.initialPressure,'X',obj.initialMoleFractions);
end
 
y0 = [intEnergy_mass(obj.gas)
    1.0/density(obj.gas)
    massFractions(obj.gas)]; %#ok<*CPROP,*PROP>
% replaces the default value of 1/10 * (t2-t1)


% define the cases
if isempty(obj.heatflux) && isempty(obj.vdot) && isempty(obj.area)
    rCases = 1; % adiabatic const volume reactor
elseif isempty(obj.heatflux) && ~isempty(obj.vdot) && ~isempty(obj.area)
    rCases = 2; % adiabatic volume changing reactor, area must also be given
elseif ~isempty(obj.heatflux) && isempty(obj.vdot) && isempty(obj.area)
    rCases = 3; % non adiabatic const volume reactor
elseif ~isempty(obj.heatflux) && ~isempty(obj.vdot) && ~isempty(obj.area)
    rCases = 4; % all "reactor functions" are given
else
    rCases = 5; % something like "undefined", eg. vdot but no area given
end
 
options = odeset('RelTol',1.e-5,'AbsTol',1.e-12,'Stats','off');
%     'BDF','on','MaxOrder',2);  
% all cases...
switch rCases
    case 1 
        % use the internal definitions for a adiabatic constant volume reactor         
        try
            fprintf('Start ode15s 0%% .')
            out = ode15s(@reactor_ode,obj.simulationTime,y0,options,obj.gas,@vdot,@area,@heatflux);  
            fprintf('.. 100%%\n')
        catch ME
           throw(ME) 
        end
    case 2 
        % reactor with heatflux, 
        % for all rCases >1 we reduce the minimal stepsize for
        % better detecting "jumping" heat fluxes etc.
               
        try
            fprintf('Start ode15s 0%% .')            
            out = ode15s(@reactor_ode,obj.simulationTime,y0,options,obj.gas,obj.vdot,obj.area,@heatflux);
            fprintf('.. 100%%\n')
        catch ME
           throw(ME) 
        end
    case 3 
        % reactor with changing volume and given area function    
        try
            fprintf('Start ode15s 0%% .')
            out = ode15s(@reactor_ode,obj.simulationTime,y0,options,obj.gas,@vdot,@area,obj.heatflux);
            fprintf('.. 100%%\n')
        catch ME
           throw(ME) 
        end
    case 4
        % reactor with the full set of optional functions     
       try
            fprintf('Start ode15s 0%% .')
            out = ode15s(@reactor_ode,obj.simulationTime,y0,options,obj.gas,obj.vdot,obj.area,obj.heatflux);
            fprintf('.. 100%%\n')
        catch ME
           throw(ME) 
        end
    case 5
        % only vdot given        
        try       
            fprintf('Start ode15s 0%% .')      
            out = ode15s(@reactor_ode,obj.simulationTime,y0,options,obj.gas,obj.vdot,@area,@heatflux);            
            fprintf('.. 100%%\n')
        catch ME          
           throw(ME)
        end
    otherwise
        fprint('You should never see this message, from now some very stange things may happen\n')
end


% simple postprocessing

obj.times = out.x ;
obj.vol = out.y(2,:);
obj.massFractions = out.y(3:end,:);
if nargin==1
    marker = 'OH';   
    fprintf(['solveForIgnition: Use species >OH< as marker for ignition.\n ']);
else
    marker = varargin{1}; 
    [~,in] = obj.lookForSpecies(marker);
    if isempty(in)
        ME = MException('solveForIgnition:MarkerNotFound',...
            ['The marker species >',marker,'< is not in the mechanism.']);
        throw(ME)
    else
        fprintf(['solveForIgnition: Use species ',marker,...
        ' as marker for ignition.\n ']);
    end
end

try 
    marker_fractions = obj.getMassFractionsByName (marker);
catch ME
    ME.identifier
    throw(ME)
end

dOH = dot(obj.times,marker_fractions);
[~,indx] = max(dOH(1,:),[],2);
obj.ignt = obj.times(1,indx(1));
% 

 


 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% some local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

% local functions, defining the reactor
    function v = vdot(~, ~, ~)
    %prameter vdot(t, vol, gas)
        v = 0.0;                                    % uncomment for constant volume
        % v = 1.e11 * (pressure(obj.gas) - 101325.0);    % holds pressure very
        % close to 1 atm
    end
%  heat flux (W/m^2).
    function q = heatflux(~,~)
        q = 0.0;                                    % adiabatic
    end
% surface area (m^2). Used only to compute heat transfer.
    function a = area(~,~)
        a = 1.0;
    end
% the first derivative of d value /d time, computed by central
% differences
    function dt = dot(time,value)
        delta_t = time(2:end)-time(1:end-1);
        dt = (value(:,2)-value(:,1))./delta_t(1);
        dt = [dt,(value(:,3:end)-value(:,1:end-2))*diag(1./(delta_t(1:end-1)+delta_t(2:end)))];
        dt = [dt,(value(:,end)-value(:,end-1))./delta_t(end)];
    end
end

