
function flag =  test_piston % like simulation
% main program for testing  a system with compression (dep. on time)
global v0
try
    br = batchReactor('/home/uwe/Arbeit/mechanismen/gri30.cti');
 
    br.initialPressure = 0.88*oneatm;
   
  
    br.initialTemperature = 490.93; 
    br.simulationTime = [0.0,1.0];
    br.initialMoleFractions = 'CH4:1,O2:3,N2:9';
    
    br.vdot = @piston;
    br.area = @area;
    br.heatflux = @heatf;
    
    set(br.gas,'T',br.initialTemperature,'P',br.initialPressure,'X',br.initialMoleFractions);    
    v0 = 1/density(br.gas);

    br.reactor = 'commonNetReactor' ;    
   
%%
    tic
    br.solve;
    toc
   
 
    figure(1)
    br.plotMassFractions(br.getLivingSpeciesNames(1e-3))
 

  
   
    flag = 1;   
catch ME
    flag = 0;
    ME.rethrow   
    throwAsCaller(ME)
end

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