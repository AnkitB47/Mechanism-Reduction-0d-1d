function dydt = conhp(obj,~,y)
% conhp method gives the reaction rate for const pressure systems.
% It is private, so U can use it only within methods of class batchReactor.
%    w = obj.conhp(time,y), 
% where time is a dummy and will be not used by conhp.
%
%    Evaluates the system of ordinary differential
%    equations for an adiabatic, constant-pressure,
%    zero-dimensional reactor. It assumes that the 'gas' object
%    represents a reacting ideal gas mixture. 
% conhp Method is inspired by conhp from CANTERA 
% (c) 2012. U.Pruefert VIRTUHCON/TUBAF

% Set the state of the gas, based on the current solution vector.

set(obj.gas, 'T', y(1), 'P', pressure(obj.gas), 'Y', y(2:end));
mw = molecularWeights(obj.gas);
% energy equation
wdot = netProdRates(obj.gas);
tdot = - temperature(obj.gas) * gasconstant * enthalpies_RT(obj.gas)' ...
       * wdot / (density(obj.gas)*cp_mass(obj.gas));
 
% species equations
rrho = 1.0/density(obj.gas);
dydt = [tdot
        rrho*mw.*wdot];
end

