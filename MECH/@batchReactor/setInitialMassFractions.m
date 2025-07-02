function obj = setInitialMassFractions(obj,names,values)
% obj = obj.setMassFractions(names,values)
% this function sets the initial mass fraction property,
% Example:
%   obj = obj.setMassFractions({'H2 O2'},[1 2]);
obj.initialMassFractions = vector2CanteraString(names,values); 
end