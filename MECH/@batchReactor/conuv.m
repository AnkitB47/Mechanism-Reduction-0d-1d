function dydt = conuv(obj,~,y)
% CONUV ODE system for a constant-volume, adiabatic reactor.
%
%    Function CONUV evaluates the system of ordinary differential
%    equations for an adiabatic, constant-volume,
%    zero-dimensional reactor. It assumes that the 'obj.gas' object
%    represents a reacting ideal obj.gas mixture. 


% Set the state of the obj.gas, based on the current solution vector.
set(obj.gas, 'T', y(1), 'Rho', density(obj.gas), 'Y', y(2:end));
nsp = nSpecies(obj.gas);
mw = molecularWeights(obj.gas);
% energy equation
wdot = netProdRates(obj.gas);
tdot = - temperature(obj.gas) * gasconstant * (enthalpies_RT(obj.gas) - ones(nsp,1))' ...
       * wdot / (density(obj.gas)*cv_mass(obj.gas));

rrho = 1.0/density(obj.gas);
% set up column vector for dydt
dydt = [ tdot 
        rrho*mw.*wdot];



end
