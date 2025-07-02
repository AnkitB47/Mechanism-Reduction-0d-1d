function fractionsOut= mixFractions(fractions1,fractions2,a)
% This function mix the fractions given by the firts two arguments,
% asssumed as a valid fractions strings.
% The third argument is the  ratio phi. We put 1/phi of the mixture in
% the first argument to the initial mass fractions in
% obj.initialMassFraction
% (C) 2011 by Uwe Pruefert for VIRTUHCON
if a>1 || a <0
    error('matlab:batchReactor:friends:mixFractions',...
        'Wrong value for mixing ratio: must be between 0 and 1.')
end
speciesListOne = getSpeciesFromInitialFractions(fractions1);
fractionsListOne = a * getMoleFractionsFromInitialFractions(fractions1)/...
    sum(getMoleFractionsFromInitialFractions(fractions1));

speciesListTwo = getSpeciesFromInitialFractions(fractions2);
fractionsListTwo =  (1-a) * getMoleFractionsFromInitialFractions(fractions2)/...
    sum(getMoleFractionsFromInitialFractions(fractions2));

% mix the species and the fractions

newSpecies = logical(1:length(speciesListTwo));

for k = 1:length(speciesListOne)    
    for l = 1:length(speciesListTwo)        
        if strcmp(speciesListOne{k},speciesListTwo{l})            
            fractionsListOne(k) =  fractionsListOne(k) + fractionsListTwo(l);
            newSpecies(l) = false;           
            break;        
        end      
    end
end

% new species is in the end a logical vector where true is when the species
% name is not fpount in the initial mass fractions of batchReactor object

speciesListOne = [speciesListOne speciesListTwo(newSpecies)];
fractionsListOne = [fractionsListOne fractionsListTwo(newSpecies)];

fractionsOut  = vector2CanteraString(speciesListOne,fractionsListOne);


end