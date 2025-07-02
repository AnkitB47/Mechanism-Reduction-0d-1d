function obj = commonReactor(obj)
obj.reactor = 'commonReactor';
obj.clearFractions;
if strcmp(obj.equilibrium,'true')
    obj. simulationTime = [0, Inf];
end


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
if ~isempty(obj.heatflux) && ~isempty(obj.vdot) && ~isempty(obj.area)
    rCases = 4; % all "reactor functions" are given
else
   ME = MException('batchReactor:commonReactor:noBC',...
       'For commonReactor heatflux, vdot and area properties must be given');
   throwAsCaller(ME)
end

switch obj.equilibrium
    case 'true'
        options = odeset('RelTol',1.e-5,'AbsTol',1.e-12,...
            'Stats','off','Events',@eventfun);
    case 'false'
        options = odeset('RelTol',1.e-5,'AbsTol',...
            1.e-12,'Stats','off');
    otherwise
end

% all cases...
switch rCases    
    case 4
        % reactor with the full set of optional functions     
       try
            fprintf('Start ode15s 0%% .')
            out = ode15s(@obj.creactor,obj.simulationTime,y0,options);           
            fprintf('.. 100%%\n')
        catch ME
            ME.message  
       end    
    otherwise
        ME = MException('batchReactor:commonReactor:noBC',...
            'For commonReactor heatflux, vdot and area properties must be given');
        throwAsCaller(ME);
end
pv = output(out,obj.gas);
dT = dot(pv(1,:),pv(2:end,:));

[~,indx] = max(dT(1,:),[],2);
obj.times = pv(1,:);
if strcmp(obj.equilibrium,'true')
    obj. simulationTime = [0, obj.times(end)];
end
obj.temp = pv(2,:);
obj.dens = pv(3,:);
obj.pres =  pv(4,:);
% 

obj.setMassFractions(pv(5:nSpecies(obj.gas)+4,:));
obj.setMoleFractions(pv(nSpecies(obj.gas)+5:end,:));

% we define the ignition delay time as the (first) maximum of dT/dt
obj.ignt = pv(1,indx(1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This functions in only for making the code cleaner. It extracs the
% different parts from the output structure of ode15s at every time step.
% Format    time 
%           temperature
%           density
%           pressure
%           massfraction species_1
%                   .
%                   :
%           massfraction species_nsp


    function [value,isterminal,direction] = eventfun(t,y)
        % event function that stopps the evaluation of the time
        % integration. Criteria for being stationary is || dY/dt(t,Y)|| <
        % 1e-3.
        value =  norm(obj.creactor(t,y),Inf)-1e-3;%  
        isterminal = 1;   % stop the integration
        direction =  -1;   % negative direction
    end
    function pv = output(s,gas)                
        soln = s.y;
        [~,n] = size(s.x);
        nSpec = nSpecies(gas);
        pv = zeros(2*nSpec + 4, n);                
        %                 set(obj.gas,'T',1001.0,'P',oneatm);
        for j = 1:n
            ss = soln(:,j);
            y = ss(3:end);  % skipp total internal 
            % energy and volumen component
            mass = sum(y);
            u_mass = ss(1)/mass; % intEnergy
            v_mass = ss(2)/mass; % 1/density
            setMassFractions(gas, y);
            setState_UV(obj.gas, [u_mass v_mass]);                    
            pv(1,j) = s.x(j);
            pv(2,j) = temperature(gas);
            pv(3,j) = density(gas);
            pv(4,j) = pressure(gas);
            pv(5:nSpec+4,j) = y; 
            pv(nSpec+5:end,j) = moleFraction(gas,speciesNames(gas))';                   
        end
    end

% local functions, defining the reactor
    function v = vdot(~, ~, ~)
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

