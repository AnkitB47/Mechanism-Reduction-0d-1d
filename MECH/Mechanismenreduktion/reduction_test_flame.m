
baseMech = '/home/uwe/data/pox-10-1.cti';


% the associated individual object
base_individual = individual(baseMech); 
%
%-----------------------------------


%% compute the referenc solution
sol = flame1d(baseMech,'gas_1'); 
sol.initialPressure = oneatm;
sol.initialTemperature = 310;
sol.mdotLeft = .3;
sol.tInLeft = 450;
sol.mdotRight = .6;
sol.tInRight = 450;

sol.compLeft = 'CH4:1';
sol.compRight = 'O2:.47,H2O:.53';
sol = sol.setSpace(linspace(0,0.02,10));

sol.flameType = 'counterFlowFlame';
sol.solve;

figure(1)
sol.plotMassFractions;
drawnow
figure(2)
plot(sol.temp);drawnow
return

%% setup the options for the GA

options = optionsGA();
options.ngen = 15;

options.selectionStrategie = 'survivors';
options.selectionThreshold = 4;
options.prob = 0.3;
options.nind = 20;


ind = individual('/home/uwe/data/gri30.cti');
 
ind = ind.fix(sol.mainSpecies);
% ind = ind.removeSpecies(1:sol.nSpecies);
%%





%% 
g = ga(ind,@dist1,sol);
g.guessFixSpecies;
g.aimsToReduce = 'species';

g.reference = sol;
g.options = options;
g.genotype.fix(ind);

g.genotype = g.genotype.removeSpecies(1:sol.nSpecies);
 
%%
g.solve(options);
%%

marker = {'--', ':' ,'-.', '-'} ;
for k = 1:3
    g.currentGeneration(k)
    sol.setMechanism( g.currentGeneration(k).mechanism);
    sol.clearFractions;
    sol.solve  
    sol.plotMassFractions([],marker{k})
    drawnow
    hold on
end


return
removeSpeciesCTI('/home/pruefert/mechanismsCO/mech85818049.cti','mech20CO.cti',sol.getDeadSpeciesIndex,sol.deadSpeciesNames)
%% 


    
        



 
