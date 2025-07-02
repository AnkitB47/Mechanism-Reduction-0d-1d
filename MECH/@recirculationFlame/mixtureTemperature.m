function T = mixtureTemperature(obj,mixture1,mixture2,T1,T2,phi)
%  T = obj.mixTemperature(mixture1,mixture2,T1,T2,phi)
% computes the temperature of a mixture with ratio phi
% mixturex must be a cantera mixture fractions string

% set up a gas for computing cp twice, for both mixtures
set(obj.reactor.gas,'Y',mixture1);
cp1 = cp_mass(obj.reactor.gas);
set(obj.reactor.gas,'Y',mixture2);
cp2 = cp_mass(obj.reactor.gas);

T = (phi*cp1*T1 + (1- phi)*cp2*T2)/(phi*cp1 + (1-phi)*cp2);

end