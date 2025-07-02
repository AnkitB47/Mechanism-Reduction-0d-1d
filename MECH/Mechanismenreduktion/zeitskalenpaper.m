function zeitskalenpaper(problem)%% clear the workspace 
mech = '/home/pruefert/data/gri30.cti';
reactor = 'constPressReactor';
switch problem
    case 1
        %  reactor 1
        mix = mixture;
        initialTemp = 1736;
        progress =  {'H2O' 'CO2' 'CO' 'H2' 'H' 'O' 'OH'}; 
    case 2
        % reactor 2. Note: These are Mole Fractions
        mix = 'CH4:1,CO2:1,O2:2';
        initialTemp = 2785;
        progress = {'H2O' 'CO2' 'CO' 'H2' 'H' 'O' 'OH'};          
    case 3        
        % reactor 3
        mix = mixture3;
        initialTemp = 2561;    
        progress = {'H2O' 'CO2' 'CO' 'H2' 'H' 'O' 'OH'};    
    case 4
        mix =  'H2:2.0,O2:1';
        initialTemp = 800;
        progress = {'H2O'};  
        mech = '/home/pruefert/data/H2_onestep.cti';
        reactor = 'constPressReactor';
    otherwise
end
 
timeToSolve = [0,3];
p0 = 1*oneatm;


% mech = '/home/pruefert/data/water_shift_9_species.cti';  
%% 
% one
br = batchReactor(mech);
br.equilibrium = 'false'; 
br.initialPressure = p0;
br.initialTemperature = initialTemp; 
br.simulationTime = timeToSolve;
if problem==2 || problem==4
    br.initialMoleFractions = mix;
else
    br.initialMassFractions = mix;
end
br.reactor = reactor ; 
 
br.solve;  
 

%%
 
%  
figure(1) 
clf
hold on
br.plotMassFractions();
ylabel('Y  ','FontSize',20)
xlabel('t [s]','FontSize',20) 
drawnow
 
%%  

 
% clf  
ts = {'spts' 'evts' 'pvts' 'ijts'};
cl = {'ko' 'rs' 'bd' 'm.'};
figure(2)
clf
for k = 1:4
    hold on
    if k == 3        
        br.pvts(progress)        
    elseif k == 2 && problem==4
        br.evts('eps',0);
    elseif k == 4 && problem==4
        br.ijts('eps',0);
    else
        br.timeScale(ts{k});
    end
    figure(2)
    loglog(br.times,min(br.lambda,[],1),cl{k});
    ylabel('\tau [s]','FontSize',20)
    xlabel('t [s]','FontSize',20)  
    %%   
    
    drawnow
    hold on
end
% 
figure(2)
legend(upper(ts),'FontSize',20); 
ylabel('\tau [s]','FontSize',20)
xlabel('t [s]','FontSize',20)


end
