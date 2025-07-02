function speciesList = getInitialSpecies(obj)
% extracts species names from the (user input) inital mole fractions            
if ~isempty(obj.initialMoleFractions)
    speciesList = ...
        obj.getSpeciesFromInitialFractions(obj.initialMoleFractions);
    return
end
if ~isempty(obj.initialMassFractions)
    speciesList = ...
        obj.getSpeciesFromInitialFractions(obj.initialMassFractions);
    return
end
error('batchReactor:getInitialSpecies:emptyList',...
    'The species list is empty.')
end 