function test_POD_Batch_Reactor
%% clear the workspace 
clc
clear all

epsSV = 1e-10;
maxSV = 30;

%% setting up the problem
 
mix = {'CH4:1, O2:4, CO2:1' ...
    'CH4:1, O2:2, CO2:1'...
   'CH4:1, O2:1, CO2:1'...
   'CH4:1, O2:4, CO:1'...
   'CH4:1, O2:2, CO:2'...
   'CH4:1, O2:1, CO:1'... 
   };
initialTemp = [850:50:1000];
timeToSolve = [0,5];
p0 = [4:2:10]*oneatm;
reactor = 'constPressReactor';
testMech = '/home/pruefert/data/gri30.cti';

%% set up two reactor objects
% one
br = batchReactor(testMech);
br.reactor = reactor ;
for k = 3:3
    br.initialPressure = p0(1);
    br.initialTemperature = initialTemp(1); 
    br.simulationTime = timeToSolve;
    br.initialMoleFractions = mix{k};
    br.solve;  
end
[V,S,~] = svd([br.temp;br.moleFractions],'econ');
idx = diag(S)>epsSV;
in = 1:length(br.times);
indx = min(maxSV,length(in(idx)));
V = V(:,indx);



figure(1)
loglog(diag(S)) 
drawnow
    function dt = dw(t,y)
        dt = V'*w(t,V*y);
    end
%% solve all objects and compute the time scales 
    function  dydt  = w(~,y)
        % source term using CANTERA
        
        
        set(gas, 'T', y(1), 'P', pressure(gas), 'Y', y(2:end));
        mw = molecularWeights(gas);
        % energy equation
        wdot = netProdRates(gas);
        tdot = - temperature(gas) * gasconstant * enthalpies_RT(gas)' ...
            * wdot / (density(gas)*cp_mass(gas));
        
        % species equations
        rrho = 1.0/density(gas);
        dydt = [tdot
            rrho*mw.*wdot];
    end


%  
%  br.plotMoleFractions
 
end
 


