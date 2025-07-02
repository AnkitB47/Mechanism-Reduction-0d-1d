% a constant volume reactor using conuv and ode15s 
function constVolumeReactor(obj)

obj.reactor = 'constVolumeReactor';
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
y0 = [temperature(obj.gas)
    massFractions(obj.gas)];
%                 ones(length(y0),length(y0))
% replace the default value of 1/10 * (t2-t1)

switch obj.equilibrium
    case 'true'
        options = odeset('RelTol',1.e-5,'AbsTol',1.e-12,...
            'Stats','off','Events',@eventfun);
    case 'false'
        options = odeset('RelTol',1.e-5,'AbsTol',1.e-12,'Stats','off');
    otherwise
end

try
    fprintf('Start ode15s 0%% .')
    out = ode15s(@obj.conuv,obj.simulationTime,y0,options);
    fprintf('.. 100%%\n')
catch ME
    ME.message  
end
pv = output(out);
dT = dot(pv(1,:),pv(2:end,:));

[~,indx] = max(dT(1,:),[],2);
obj.times = pv(1,:);
if strcmp(obj.equilibrium,'true')
    obj. simulationTime = [0, obj.times(end)];
end
obj.temp = pv(2,:);
obj.dens = pv(3,:);
obj.pres =  pv(4,:); 
obj.vol = 1./obj.dens;
% 
obj.setMassFractions(pv(5:nSpecies(obj.gas)+4,:));
obj.setMoleFractions(pv(nSpecies(obj.gas)+5:end,:));
% we define the ignition delay time as the (first) maximum of dT/dt
obj.ignt = pv(1,indx(1));


    function [value,isterminal,direction] = eventfun(t,y)
        % event function that stopps the evaluation of the time
        % integration. Criteria for being stationary is || dY/dt(t,Y)|| <
        % 1e-3.
        value =  norm(obj.conuv(t,y),Inf)-1e-3;%  
        isterminal = 1;   % stop the integration
        direction =  -1;   % negative direction
    end
% output
    function pv = output(s)                
        soln = s.y;
        [~,n] = size(s.x);                
        nSpec = nSpecies(obj.gas);
        pv = zeros(2*nSpec+4, n);                
        %                 set(obj.gas,'T',1001.0,'P',oneatm);
        for j = 1:n
            y = soln(2:end,j);
            mass = sum(y);
            u_mass = intEnergy_mass(obj.gas)/mass;
            v_mass = 1.0/density(obj.gas)/mass;
            setMassFractions(obj.gas, y);
            setState_UV(obj.gas, [u_mass v_mass]); 
            
            pv(1,j) = s.x(j); % time
            pv(2,j) = s.y(1,j); % temperature
            pv(3,j) = density(obj.gas);
            pv(4,j) = pressure(obj.gas);
            
            pv(5:nSpec+4,j) = y./mass; 
            pv(nSpec+5:end,j) = moleFraction(obj.gas,speciesNames(obj.gas))';      
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

