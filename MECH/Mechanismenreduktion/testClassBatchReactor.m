%% clear the workspace 
clear all
%% setting up the problem
%  mix= 'CH4:0.3193,O2:0.4246,H2O:0.2561';
% mix2 = 'CH4:0.0284,O2:0.2264,N2:0.7452';

% reactor 1
% mix = mixture;
% initialTemp = 1736;

% reactor 3
% mix = mixture3;
% initialTemp = 2561;


% reactor 2. Note: These are Mole Fractions
% mix = 'CH4:1,CO2:1,O2:2';
% initialTemp = 2785;

% mix = 'CH4:0.5,O2:1';
% initialTemp = 886;
mix = 'CO:0.2,CO2:0.3,H2:0.0,H2O:0.3';
initialTemp = 1210;

timeToSolve = [0,1000];
p0 = 1*oneatm;
reactor = 'constVolumeNetReactor';
% mech = '/home/pruefert/data/gri30.cti';
mech = 'gri30.cti';
%% 
% one
br = batchReactor(mech);
br.equilibrium = 'false'; 
br.initialPressure = p0;
br.initialTemperature = initialTemp; 
br.simulationTime = timeToSolve;
br.initialMassFractions = mix;
br.reactor = reactor ; 
% m = br.equilibrate('UV');
% 
% figure(1)
% bar(m);
br.solve;  
% figure(2)
% br.plotMassFractions;
% figure(3)
% bar(br.massFractions(:,end))

%%
figure(1)
% clf
br.plotTemp
drawnow
hold on
%  
figure(2) 
% clf
hold on
br.plotMassFractions([],'k-');
drawnow
return
%%  

hold on
% clf  
ts = {'spts' 'evts' 'pvts' 'ijts'};
cl = {'ko' 'rs' 'bd' 'm.'};
for k = 1:4
%     figure(k)
hold on
    br.timeScale(ts{k});
    loglog(br.times,min(br.lambda,[],1),cl{k});
    ylabel('\tau [s]','FontSize',20)
    xlabel('t [s]','FontSize',20)  
    %%   
    
    drawnow
    hold on
end
% 
% legend(upper(ts),'FontSize',20); 
% ylabel('\tau [s]','FontSize',20)
% xlabel('t [s]','FontSize',20)
