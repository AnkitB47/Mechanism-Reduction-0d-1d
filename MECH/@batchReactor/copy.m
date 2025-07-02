function obj2 = copy(obj1)
% copy function for batchReactor
% Because of chemistry is an handle, we need also a hardcopy method...
%  

if ~obj1.isvalid
    ME = MException('batchReactor:copy',...
        'Invalid or deleted object');
    throw(ME)
end

% Create a NEW object... 
obj2 = batchReactor();

% Copy all properties into the new object...
obj2.gas = obj1.gas;
obj2.mechanism = obj1.mechanism ;
obj2.speciesNames = obj1.speciesNames;
obj2.times = obj1.times;
obj2.lambda = obj1.lambda;
obj2.ignt = obj1.ignt;
obj2.temp = obj1.temp;
obj2.dens = obj1.dens;
obj2.pres = obj1.pres;
obj2.vol = obj1.vol;
obj2.reactor = obj1.reactor;
obj2.initialMoleFractions = obj1.initialMoleFractions;
obj2.initialMassFractions = obj1.initialMassFractions;
obj2.initialTemperature = obj1.initialTemperature;
obj2.initialPressure = obj1.initialPressure;
obj2.simulationTime = obj1.simulationTime;
obj2.equilibrium = obj1.equilibrium;
obj2.heatflux = obj1.heatflux;
obj2.vdot = obj1.vdot;
obj2.area = obj1.area;
obj2.massFractions = obj1.massFractions;
obj2.moleFractions = obj1.moleFractions;
end

