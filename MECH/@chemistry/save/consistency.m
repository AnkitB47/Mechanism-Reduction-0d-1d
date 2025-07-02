function [list,indx] = consistency(obj,list,indx)
% [list,indx] = consistency(obj,list,indx)
% checks if a species from the list is in the inital mole fractions. These
% species can not removed from the mechanism and will be removed from the
% list

speciesList = selectSpeciesFromList(obj.initalMoleFractions);





function speciesList = selectSpeciesFromList(str) 
    % extracts species names from the (user input) inital mole fractions
    % here we use it for selecting the three body from their definiton
    % sequence what is quite similar to the cantera systax of the mole
    % fractions
    speciesList = [];
    k = 1; % zaehler fuer die position
    in = true;
    species = [];
    for i=1:length(str)        
        if strcmp(str(i),':')
            in = false; 
            speciesList{k} = species;
            species = [];
            k = k+1; % next  
        elseif isLetter(str(i))
            in = true;
        end
        if allowedChar(str(i)) && in
            species = [species,str(i)];             
        end                         
    end
    end

end