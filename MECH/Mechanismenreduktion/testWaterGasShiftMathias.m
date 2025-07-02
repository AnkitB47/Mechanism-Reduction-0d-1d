%% clear the workspace 


for  initialTemp = 1200
mix = 'CO:0.2,CO2:0.3,H2:0.0,H2O:0.3'; 

timeToSolve = [0,10000];
p0 = 1*oneatm;
reactor = 'constPressNetReactor';

%% 

% reference
mech = '/home/uwe/data/gri30.cti';
ref = batchReactor(mech);
ref.equilibrium = 'true'; 
ref.initialPressure = p0;
ref.initialTemperature = initialTemp; 
ref.simulationTime = timeToSolve;
ref.initialMassFractions = mix;
ref.reactor = reactor ; 
ref.solve; 
% reduced mech
mech = '/home/uwe/data/water_shift_9_species.cti'; 
br = batchReactor(mech);
br.equilibrium = 'false'; 
br.initialPressure = p0;
br.initialTemperature = initialTemp; 
br.simulationTime = ref.simulationTime;
br.initialMassFractions = mix;
br.reactor = reactor ; 
br.solve;  

 
figure(1)
clf
br.plotTemp
drawnow
hold on
ref.plotTemp
%  
figure(2) 
clf
hold on
br.plotMassFractions(br.speciesNames);
ref.plotMassFractions(br.speciesNames);

drawnow
br.dist1(ref)
end
