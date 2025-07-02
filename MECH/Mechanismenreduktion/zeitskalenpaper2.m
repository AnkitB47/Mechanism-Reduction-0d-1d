function zeitskalenpaper2(problem)
    % clear the workspace 
    clc
    progress{1} =  {'CO2' 'H2O'};      
    progress{2} =  {'CO2' 'H2O' 'CO' 'H2'};     
    progress{3} =  {'CO2' 'H2O' 'CO' 'H2' 'OH' 'O' 'H'}; 
     
   
    timeToSolve = [0,3];
    p0 = 1*oneatm;
    reactor = 'constTempReactor';
    mech = '/home/pruefert/data/gri30.cti';
    br = batchReactor(mech);
    br.equilibrium = 'false'; 
    br.initialPressure = p0;
    
  
    
    
    br.simulationTime = timeToSolve;    
    
    switch problem
        case 1
            %  reactor 1
            mix = mixture;
            initialTemp = 1736;
            br.initialMassFractions = mix;
        case 2
            % reactor 2. Note: These are Mole Fractions
            mix = 'CH4:1,CO2:1,O2:2';
            initialTemp = 2785;
            br.initialMoleFractions = mix;
        case 3        
            % reactor 3
            mix = mixture3;
            initialTemp = 2561;    
            br.initialMassFractions = mix; 
        otherwise
        %
    end
   
    br.initialTemperature = initialTemp;  
    br.reactor = reactor ;
    
    br.solve;  
    
    %  
    figure(1) 
    clf
    hold on
    br.plotMassFractions();
    ylabel('Y  ','FontSize',20)
    xlabel('t [s]','FontSize',20) 
    drawnow 
    figure(2)
    clf
    cl = {'ko' 'rs' 'bd' 'm.'}; 
    for k = 1:3
        hold on
        
        br.pvts(progress{k})        
        
        
        figure(2)
        loglog(br.times,min(br.lambda,[],1),cl{k});
        ylabel('\tau [s]','FontSize',20)
        xlabel('t [s]','FontSize',20)  
        %    
        
        drawnow
        hold on
    end
    br.spts;
    figure(2)
    loglog(br.times,min(br.lambda,[],1),cl{4});
    ylabel('\tau [s]','FontSize',20)
    xlabel('t [s]','FontSize',20)  
    %    
    
    drawnow
    hold on
    % 
    figure(2)
    legend({'PV1' 'PV2' 'PV3' 'SPTS'},'FontSize',20); 
    ylabel('\tau [s]','FontSize',20)
    xlabel('t [s]','FontSize',20)
end