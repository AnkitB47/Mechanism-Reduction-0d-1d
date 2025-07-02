function obj = setMoleFractions(obj,names,values)
% obj = obj.importFractions(names,values)
obj.initialMoleFractions = vector2CanteraString(names,values); 
end
