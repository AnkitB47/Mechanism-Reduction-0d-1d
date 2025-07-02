function obj = constVolumeNetReactor(obj)

% reset and mark the object
obj.reactor = 'constVolumeNetReactor';
obj.clearFractions;

% set the initial conditions

if ~isempty(obj.initialMoleFractions)
    set(obj.gas,'T',obj.initialTemperature,'P',obj.initialPressure,'X',obj.initialMoleFractions);
elseif ~isempty(obj.initialMassFractions)
    set(obj.gas,'T',obj.initialTemperature,'P',obj.initialPressure,'Y',obj.initialMassFractions);
else
    error('batchReactor:constVolumeNetReactor:noFractionsGiven',...
        'The reactor must be initialized eigher with mole or mass frations.')
end
% create a reactor, and insert the obj.gas
r = Reactor(obj.gas);

% create a reactor network and insert the reactor
network = ReactorNet({r});

% N steps: for very stiff problems, 800 seems to be sufficient, but for
% non-stiff problems it can be a little bit overkill.
N = 800;

% simulation time, not neccessary [0,t]...
t = min(obj.simulationTime);
dt = (max(obj.simulationTime)-min(obj.simulationTime))/N;

% clear all fields in the object 
obj.clearFractions;
obj.times = [];         obj.temp = [];
obj.dens = [];          obj.pres = [];

fprintf('Solve reactor net...')  
for n = 1:N
    t = t + dt;
    try
        advance(network, t);
    catch ME
        % hmmmm????
        ME.rethrow;
    end
    obj.times(n) = time(network);
    
    obj.temp(n) = temperature(r);
    obj.dens(n) = density(r);
    obj.pres(n) = pressure(r);
    
    obj.moleFractions(:,n) = moleFractions(obj.gas);
    obj.massFractions(:,n) = massFractions(obj.gas);       
end
fprintf(' done\n')

 
% ignition time as argmax(dT/dt)
[~,indx] = max((obj.temp(2:end)-obj.temp(1:end-1))./(obj.times(2:end)-obj.times(1:end-1)));
obj.ignt = obj.times(indx);
end