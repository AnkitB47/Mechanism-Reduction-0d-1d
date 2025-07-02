function  constTempReactor(obj)

obj.reactor = 'constTempReactor';
obj.clearFractions;
if strcmp(obj.equilibrium,'true')
    obj. simulationTime = [0, Inf];
end


if  isa(obj.mechanism,'char') && isempty(obj.gas)
    obj.gas = IdealGasMix(obj.mechanism)  ;
    fprintf('chemistry.constTempReactor: Make obj.gas object\n')
elseif ~isIdealGas(obj.gas)   
    error('obj.gas is not an Ideal Gas Mix')
end


if isempty(obj.initialMoleFractions)
    set(obj.gas,'T',obj.initialTemperature,'P',obj.initialPressure);                
    set(obj.gas,'Y',obj.initialMassFractions);                
else
    set(obj.gas,'T',obj.initialTemperature,'P',obj.initialPressure,'X',obj.initialMoleFractions);
end

switch obj.equilibrium
    case 'true'
        options = odeset('RelTol',1.e-5,'AbsTol',1.e-12,...
            'Stats','off','Events',@eventfun);
    case 'false'
        options = odeset('RelTol',1.e-5,'AbsTol',1.e-12,'Stats','off');
    otherwise
end

% here only fractions, there is no temperature in the solution
y0 =  massFractions(obj.gas);
try
    fprintf('Start ode15s 0%% .')
    out = ode15s(@obj.conht,obj.simulationTime,y0,options);
    fprintf('.. 100%%\n')
catch ME
    ME.message      
end
pv = output(out);
dT = dot(pv(1,:),pv(2,:));

[~,indx] = max(dT(1,:));

obj.times = pv(1,:);
if strcmp(obj.equilibrium,'true')
    obj. simulationTime = [0, obj.times(end)];
end
obj.temp = temperature(obj.gas)*ones(size(obj.times));
obj.dens = pv(2,:); 
obj.vol = 1.0./obj.dens;
obj.pres = pv(3,:);
% 
 obj.setMassFractions(pv(4:nSpecies(obj.gas)+3,:));
 obj.setMoleFractions(pv(nSpecies(obj.gas)+4:end,:));
% we define the ignition delay time as the (first) maximum of dT/dt
obj.ignt = pv(1,indx(1));%             

% output

    function [value,isterminal,direction] = eventfun(t,y)
        % event function that stopps the evaluation of the time
        % integration. Criteria for being stationary is || dY/dt(t,Y)|| <
        % 1e-3.
        value =  norm(obj.conht(t,y),Inf)-1e-3;%  
        isterminal = 1;   % stop the integration
        direction =  -1;   % negative direction
    end
    function pv = output(s)                
        soln = s.y;
        [~,n] = size(s.x);         
        pv = zeros(2*obj.nSpecies + 3, n);            
       
        for j = 1:n
            y = soln(:,j);
            setMassFractions(obj.gas,y);             
            pv(1,j) = s.x(j); % time
            pv(2,j) = density(obj.gas);    
            pv(3,j) = pressure(obj.gas);
            pv(4:obj.nSpecies+3,j) = y; 
            pv(obj.nSpecies+4:end,j) = moleFraction(obj.gas,obj.speciesNames)';      
        end
    end
% derivative overwrites DOT operator locally
    function dt = dot(time,value)
        delta_t = time(2:end)-time(1:end-1);
        dt = (value(:,2)-value(:,1))./delta_t(1);
        dt = [dt,(value(:,3:end)-value(:,1:end-2))*diag(1./(delta_t(1:end-1)+delta_t(2:end)))];
        dt = [dt,(value(:,end)-value(:,end-1))./delta_t(end)];            
    end
end