function b = isDeadReaction(reactionLine,deads)
% b= isDeadReaction(reactionLine) 
%   decides if an reaction is dead, i.e. if it contains a dead species
b = false;
species = chemistry.selectSpeciesFromList(reactionLine);
for k = 1:length(deads)
    for l = 1:length(species)
        if strcmp(deads{k},species{l})
            b = true;
            return
        end
    end
end
            

end

