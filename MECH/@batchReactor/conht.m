function dydt = conht(obj,~,y)
% CONUV ODE system for a constant-volume, adiabatic reactor.
%
%    Function CONUV evaluates the system of ordinary differential
%    equations for an adiabatic, constant-volume,
%    zero-dimensional reactor. It assumes that the 'gas' object
%    represents a reacting ideal gas mixture. 


% Set the state of the gas, based on the current solution vector.
set(obj.gas,  'Rho', density(obj.gas), 'Y', y);
mw = molecularWeights(obj.gas);

% energy equation
wdot = netProdRates(obj.gas);
rrho = 1.0/density(obj.gas);
% set up column vector for dydt
dydt =  rrho*mw.*wdot;



end
