clear classes
clear all
clc
clf
fl = recirculationFlame('gri30.cti','gri30_mix');
fl.flame.compLeft = fl.reactor.mass2mole('CH4:0.278,CO2:0.458,H2O:0,H2:0.07,CO:0.194');
fl.flame.compRight = 'O2:1';
fl.flame.space = linspace(0,.05,20);
fl.flame.initialTemperature = 300;
fl.flame.initialPressure =  oneatm;
fl.flame.mdotLeft = 0.1;
fl.flame.mdotRight = 0.2;
fl.flame.tInLeft = 300;
fl.flame.tInRight = 300; 

fl.reactor.initialPressure=oneatm;
fl.reactor.simulationTime = [0,100];
fl.reactor.reactor = 'constTempReactor';
fl.recirculatingMixture = fl.flame.compLeft;
% fl.reactor.initialTemperature = 1300;
fl.ratio = 18/39;

%% 
fl = fl.solve;    

%%


%%
hold on
figure(1)

fl.plotMassFractions()
drawnow
figure(2)
hold on
fl.plotTemp;

%%
fl.reactor.speciesTimeScales(5,fl.reactor.getInitialSpecies);
semilogy(fl.reactor.times(2:end),(fl.reactor.times(2:end)-fl.reactor.times(1:end-1)),'.')
%%
fl.reactor.eigenvalueTimeScales('step',2);
semilogy(fl.reactor.times(2:end),(fl.reactor.times(2:end)-fl.reactor.times(1:end-1)),'.')


