

%%
mix = 'CO:0.2,CO2:0.3,H2:0.0,H2O:0.3';
initialTemp = 2000;

timeToSolve = [0,1];
p0 = 1*oneatm;
reactor = 'constVolumeNetReactor';
mech = '/home/uwe/data/gri30.cti';
%% 
% one
sol = batchReactor(mech);
sol.equilibrium = 'false'; 
sol.initialPressure = p0;
sol.initialTemperature = initialTemp; 
sol.simulationTime = timeToSolve;
sol.initialMassFractions = mix;
sol.reactor = reactor ; 
 

sol.solve;  
figure(1)
hold on
sol.plotTemp;
drawnow
%  
figure(2) 
% clf
hold on
sol.plotMassFractions([],'k-');
drawnow 

hold on


%%
ind = individual('/home/uwe/data/gri30.cti');
 
ind = ind.fix(sol.mainSpecies);
% ind = ind.removeSpecies(1:sol.nSpecies);
%%
options = optionsGA;
options.ngen = 4;
options.prob = 0.3;
options.nind = 20;
options.selectionThreshold = 4;
%%
% load g_mathias_save_4rd_run

% 
g = ga(ind,@dist21,sol);
g.guessFixSpecies;
g.aimsToReduce = 'species';

g.reference = sol;
g.options = options;
g.genotype.fix(ind);

g.genotype = g.genotype.removeSpecies(1:sol.nSpecies);
 

g.solve(options);
%%
figure(1);
figure(2);
for k = 1:2
    solOpt = batchReactor(g.currentGeneration(k));
    solOpt.initialPressure = 1*oneatm;
    solOpt.initialTemperature = sol.initialTemperature; 
    solOpt.simulationTime = [0 1.5];
    solOpt.initialMassFractions = 'CO:0.2,CO2:0.3,H2:0.0,H2O:0.3'  ;
    solOpt.reactor = 'constPressReactor';
    solOpt.solve;  
    sol.dist21(solOpt)
    figure(1)
    hold on
    solOpt.plotTemp();
    figure(2)
    hold on
    solOpt.plotMassFractions();
    hold on
    drawnow
end
% save('g_HB_4.mat','g')
% exit