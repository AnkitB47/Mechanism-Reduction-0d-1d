br = batchReactor('/home/pruefert/data/gri30.cti');
br.initialPressure = oneatm;
br.initialTemperature = 1752.93; 
br.simulationTime = [0,3];
br.initialMassFractions = mixture;
br.reactor = 'constTempReactor' ; 
br = br.solve;  
br =  br.speciesTimeScales('eps',1e-9);
%%
br1 = br;
br1 = br1.setMechanism ('/home/pruefert/data/grimech30_113.cti');
br1 = br1.solve;  
br1 =  br1.speciesTimeScales('eps',1e-9);
br.dist12(br1)