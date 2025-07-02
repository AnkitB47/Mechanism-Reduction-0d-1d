
function flag =  test_piston % like simulation
% main program for testing  a system with compression (dep. on time)
global v0
% try
    mech = '/usr/local/share/cantera/data/gri30.cti';
    br = batchReactor(mech);
 
    br.initialPressure = 0.88*oneatm;
   
  
    br.initialTemperature = 490.93; 
    br.simulationTime = [0.0,1.0];
    br.initialMoleFractions = 'CH4:1,O2:3,N2:9';
    
    br.vdot = @piston;
    br.area = @area;
    br.heatflux = @heatf;
    
    set(br.gas,'T',br.initialTemperature,'P',br.initialPressure,'X',br.initialMoleFractions);    
    v0 = 1/density(br.gas);

    br.reactor = 'commonReactor' ;    
   
%%
    tic
      br.solve;
    toc
%     tic
%     br = br.solveForIgnition('OH');
%     toc
%     tic
%      br.timeScaleDFS('eps',1e-3);
%     toc
 
    figure(1)
    br.plotMassFractions(br.getLivingSpeciesNames(1e-3))
 
%     figure(2)
%     br.plotPress
%     
%     figure(3)
%     br.plotTemp
%      
%    

%     figure(2) 
%     loglog(br.times,br.lambda,'b.')
%     
%     figure(4)
%     plot(br.times,br.vol)
   
%     flag = 1;   
% % catch ME
%     flag = 0;
%     warning(ME.message);   
%     throwAsCaller(ME)
% end

% local functions


    function dv = piston(t,~,~)
        %prameter vdot(t, vol, gas)
        % simulates a piston with v(0) = 1 v(pi) = 0.2
        % phi is the frequence d is the relation Vmin/Vmax
        phi = 10; 
        d = 10;       
        a = 0.5*(v0-v0/d);
        dv = -a*phi*sin(phi*t);                                   
    end

    function q = heatf(~,gas)
        q  = (400-temperature(gas))*0.3*1e3;
    end

    function ar = area(t,~)
       phi = 10; 
        d = 10;       
        a = 0.5*(v0-v0/d);           
        ar = a*(cos(phi*t)+1)+v0/d;       
    end


end