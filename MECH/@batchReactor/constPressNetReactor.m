function constPressNetReactor(obj)
% constPressNetReactor  
% Zero-dimensional kinetics: adiabatic, constant pressure.
%
%
% 
%
%    This reactor uses the class 'Reactor' for
%    zero-dimensional kinetics simulations. Here the parameters are
%    set so that the reactor is adiabatic and very close to constant
%    pressure.
%
%   (c) Uwe Pruefert for VIRTUCHON  2011


obj.reactor = 'constPressNetReactor';
% set the initial conditions
obj.clearFractions;
if ~isempty(obj.initialMoleFractions)
    set(obj.gas,'T',obj.initialTemperature,'P',obj.initialPressure,'X',obj.initialMoleFractions);
elseif ~isempty(obj.initialMassFractions)
    set(obj.gas,'T',obj.initialTemperature,'P',obj.initialPressure,'Y',obj.initialMassFractions);
else
    error('batchReactor:constPressNetReactor:noFractionsGiven',...
        'The reactor must be initialized eigher with mole or mass frations.')
end

% create a reactor, and insert the gas
r = Reactor(obj.gas);

% create a reservoir to represent the environment
set(obj.gas,'T',obj.initialTemperature,'P',obj.initialPressure,'X','N2:1');

env = Reservoir(obj.gas);

% Define a wall between the reactor and the environment and
% make it flexible, so that the pressure in the reactor is held
% at the environment pressure.
w = Wall;
install(w,r,env);

% set expansion parameter. dV/dt = KA(P_1 - P_2)
setExpansionRateCoeff(w, 1.0e6);

% set wall area
setArea(w, 1.0);

% create a reactor network and insert the reactor:
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
        advance(network,t);
    catch ME    
        ME.throw;
    end
    obj.times(n) = time(network);    
    obj.temp(n) = temperature(r);
    obj.dens(n) = density(r);
    obj.pres(n) = pressure(r);   
    obj.moleFractions(:,n) = moleFractions(obj.gas);
    obj.massFractions(:,n) = massFractions(obj.gas);
end
fprintf(' done.\n')
 
% ignition time as argmax(dT/dt)
[~,indx] = max((obj.temp(2:end)-obj.temp(1:end-1))./(obj.times(2:end)-obj.times(1:end-1)));
obj.ignt = obj.times(indx);

end
